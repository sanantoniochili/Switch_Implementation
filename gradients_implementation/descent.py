import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from cysrc.potential import *
from ase.geometry import get_distances
from linmin import bisection_linmin

import shutil
COLUMNS = shutil.get_terminal_size().columns

def prettyprint(dict_):
	import pprint
	np.set_printoptions(suppress=True)
	words = ""
	for key, value in dict_.items():
		if isinstance(value,int) or isinstance(value,float):
			words += key+" "+str(value)+" "
		else:
			print("\n", key)
			print(value)
	print(words.center(COLUMNS,"-"))


def GD(grad, **kwargs):
	"""Gradient Descent updating scheme. Returns
	the direction vector.

	"""
	return -grad/np.linalg.norm(grad)


def CG(grad, **kwargs):
	"""Conjugate Gradient updating scheme. We are
	using Polak-Ribière updating for beta factor.
	All matrices are reshaped and addressed as vectors.

	"""

	residual = -np.reshape(grad,(grad.shape[0]*grad.shape[1],1))
	last_residual = np.reshape(kwargs['residual'], 
		(kwargs['residual'].shape[0]*kwargs['residual'].shape[1],1))
	last_direction = np.reshape(kwargs['direction'], 
		(kwargs['direction'].shape[0]*kwargs['direction'].shape[1],1))
	beta = residual.T @ (residual-last_residual) / (residual.T @ residual)

	ndirection = residual + beta*last_direction
	return np.reshape(ndirection, (kwargs['direction'].shape[0],kwargs['direction'].shape[1]))


class Descent:
	def __init__(self, ftol=0.00001, gtol=0.001, 
		tol=0.00000001, iterno=1000):
		self.iterno = iterno
		self.ftol = ftol
		self.gtol = gtol
		self.tol = tol
		self.iters = 0
		self.methods = []
		self.CG_iterno = 0
		self.cattol = 0.5

	def completion_check(self, init_energy, energy, grad):
		de = abs(init_energy-energy)
		if de <= self.ftol:
			print("Iterations: {} Energy: {} Final Energy difference: {} Tolerance: {}".format(
					self.iters,energy,de,self.ftol))
			return True
		gnorm = np.linalg.norm(grad)
		if gnorm <= self.gtol:
			print("Iterations: {} Final Gnorm: {} Tolerance: {}".format(
					self.iters,gnorm,self.gtol))
			return True
		return False

	def iter_step(self, atoms, potentials, cell_wrap, last_iter={}, 
		step_func=bisection_linmin, direction_func=GD):
		"""Updating iteration step. The energy is calculated on
		new ion positions that are not inflicted to the Atoms
		object.

		Returns:
		- Direction 	: The direction vector used for the updating step
		- Positions 	: The updated ions' positions
		- Gradient 		: The gradient vector w.r.t. the positions before change
		- Energy 		: The value of energy w.r.t. the updated positions
		
		"""
		pos = np.array(atoms.positions)
		vects = np.array(atoms.get_cell())
		N = len(atoms.positions)

		ri=-last_iter['Gradient']
		pi=last_iter['Direction']

		# Gradient
		grad_coul = np.array(
			potentials['Coulomb'].calc_drv(
			pos_array=pos, vects_array=vects, N_=N))
		grad_buck = np.array(
			potentials['Buckingham'].calc_drv(
			pos_array=pos, vects_array=vects, N_=N))
		grad = grad_coul+grad_buck
		if "Lagrangian" in potentials:
			grad += np.array(
			potentials['Lagrangian'].calc_constrain_drv(
			pos_array=pos, N_=N))

		# Direction
		# If CG is used N times then change to GD for one iter
		if (direction_func.__name__=="CG"):
			if self.CG_iterno == N:
				direction_func = GD
				self.CG_iterno = 0
			else:
				self.CG_iterno += 1 # count CG iters

		pi_1 = direction_func(grad, residual=ri, direction=pi)
		
		# Add name of used method to list
		self.methods += [direction_func.__name__]

		# Calculate step size
		(step, energy) = step_func(atoms, pi_1, potentials, grad, 
			last_iter, cell_wrap, self.ftol)

		# Calculate new point on energy surface
		self.iters += 1
		# pos_temp = atoms.positions + step*pi_1
		pos_temp = cell_wrap.move_ions(atoms.positions,
			step*pi_1, vects, len(atoms.positions))

		return {'Positions':pos_temp, 'Direction':pi_1, 'Gradient':grad, 
		'Iter':self.iters, 'Step':step, 'Energy':energy}

	def repeat(self, init_energy, atoms, potentials, 
		cell_wrap, step_func=bisection_linmin, direction_func=GD):

		pos = np.array(atoms.positions)
		vects = np.array(atoms.get_cell())
		N = len(atoms.positions)

		prettyprint({'Initial Energy':init_energy})
		input()

		# First iteration (Steepest Descent) ##########################
		grad_coul = np.array(
			potentials['Coulomb'].calc_drv(
			pos_array=pos, vects_array=vects, N_=N))
		grad_buck = np.array(
			potentials['Buckingham'].calc_drv(
			pos_array=pos, vects_array=vects, N_=N))
		grad = grad_coul+grad_buck
		if "Lagrangian" in potentials:
			grad += np.array(
			potentials['Lagrangian'].calc_constrain_drv(
			pos_array=pos, N_=N))
		
		p = GD(grad)
		# Add name of used method to list
		self.methods += ['GD']
		
		# Calculate step size
		(step, energy) = step_func(atoms, p, potentials, grad,
			{'Energy': init_energy, 
			'Gradient': grad,
			'Direction': p}, cell_wrap, self.ftol)

		# Calculate new point on energy surface
		self.iters += 1
		# atoms.positions = atoms.positions + step*p
		atoms.positions = cell_wrap.move_ions(
			atoms.positions, step*p, vects, len(atoms.positions))

		# Catastrophe check
		if potentials['Buckingham'].catastrophe_check(\
			atoms.positions, self.cattol, 1) is not None:
			print("!! Detected ions too close !!".center(COLUMNS))

		iteration = {'Positions':atoms.positions, 'Direction':p, 
		'Gradient':grad, 'Iter':self.iters, 'Step':step,  'Energy':energy}

		prettyprint(iteration)
		input()

		if self.completion_check(init_energy, iteration['Energy'], iteration['Gradient']):
			return iteration

		# Rest of iterations ##########################################
		for i in range(self.iterno):
			last_iteration = iteration
			iteration = self.iter_step(atoms, potentials, 
				cell_wrap=cell_wrap,
				last_iter=last_iteration, step_func=step_func, 
				direction_func=direction_func)
			
			# Keep the newly found energy value
			self.emin = iteration['Energy']
			# Change point on PES and try a whole step
			atoms.positions = iteration['Positions']
			# Catastrophe check
			if potentials['Buckingham'].catastrophe_check(\
				atoms.positions, self.cattol, 1) is not None:
				print("!! Detected ions too close !!".center(COLUMNS))

			if self.completion_check(last_iteration['Energy'], 
				iteration['Energy'], iteration['Gradient']):
				break

			prettyprint(iteration)
			input()

		return iteration


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'-i', metavar='--input', type=str,
		help='.cif file to read')
	args = parser.parse_args()
	atoms = aread(args.i)