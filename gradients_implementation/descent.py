import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from cysrc.potential import *

def prettyprint(dict_):
	import pprint
	np.set_printoptions(suppress=True)
	for key, value in dict_.items():
		print(key)
		print(value,"\n")


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


def bisection_linmin(atoms, direction, potentials, grad, 
	init_iter, c1=0.0001, c2=0.1):

	vects = np.array(atoms.get_cell())
	N = len(atoms.positions)
	step = 1

	while(True):
		# Check if energy is dropping
		pos_temp = atoms.positions + step*direction
		energy = potentials['Coulomb'].calc(
					pos_array=pos_temp, 
					vects_array=vects, N_=N)['Electrostatic'] + \
				potentials['Buckingham'].calc(
					pos_array=pos_temp, 
					vects_array=vects, N_=N)
		if "Lagrangian" in potentials:
			energy += potentials['Lagrangian'].calc_constrain(
				pos=pos_temp, N=N)

		print("Initial energy: ",init_iter['Energy'],
			" New energy: ",energy, "Step: ",step)

		vdirection = np.reshape(direction, 
			(direction.shape[0]*direction.shape[1],1))
		vgrad = np.reshape(grad, 
			(grad.shape[0]*grad.shape[1],1))

		# Wolfe conditions
		# Make sure that the new step is large enough to decrease energy
		if (energy-init_iter['Energy']) <= step*c1*(vdirection.T @ vgrad):
			break

		# Update step size
		step = step/2
		print("New step: ",step)
	return (step, energy)


class Descent:
	def __init__(self, ftol=0.0001, gtol=0.001, 
		tol=0.00000001, iterno=1000):
		self.iterno = iterno
		self.ftol = ftol
		self.gtol = gtol
		self.tol = tol

	def iter_step(self, atoms, potentials, last_iter={}, 
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
		pi_1 = direction_func(grad, residual=ri, direction=pi)

		# Calculate step size
		(step, energy) = step_func(atoms, pi_1, potentials, grad, 
			last_iter, self.ftol)

		# Calculate new point on energy surface
		pos_temp = atoms.positions + step*pi_1

		return {'Direction':pi_1, 'Gradient':grad, 'Step':step, 
		'Positions':pos_temp, 'Energy':energy }

	def repeat(self, init_energy, atoms, potentials, 
		step_func=bisection_linmin, direction_func=GD):

		pos = np.array(atoms.positions)
		vects = np.array(atoms.get_cell())
		N = len(atoms.positions)

		# First iteration
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
		
		gnorm = np.linalg.norm(grad)
		p = -grad/gnorm
		
		# Calculate step size
		(step, energy) = step_func(atoms, p, 
			potentials, grad,
			{'Energy': init_energy, 
			'Gradient': grad,
			'Direction': p}, self.ftol)

		# Calculate new point on energy surface
		atoms.positions = atoms.positions + step*p

		iteration = {'Direction':p, 'Gradient':grad, 'Step':step, 
		'Positions':atoms.positions, 'Energy':energy}

		prettyprint(iteration)
		input()

		de = abs(init_energy-iteration['Energy'])
		if de < self.ftol:
			return iteration

		for i in range(self.iterno):
			last_iteration = iteration
			iteration = self.iter_step(atoms, potentials,
				last_iter=last_iteration, step_func=step_func, 
				direction_func=direction_func)
			
			# Keep the newly found energy value
			self.emin = iteration['Energy']
			# Change point on PES and try a whole step
			atoms.positions = iteration['Positions']

			de = abs(last_iteration['Energy']-iteration['Energy'])
			if de < self.ftol:
				print("Energy:",iteration['Energy'],"Energy difference: ",de,
					" Tolerance: ",self.ftol)
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