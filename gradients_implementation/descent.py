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
	using Polak-Ribière updating for beta factor. Beta
	array contains one scalar for each ion coordinate
	aka each optimised vector in the 2D matrix being 
	optimised.
	Returns the direction vector.

	"""

	residual = -grad
	ndirection = np.ones(grad.shape)
	beta = np.ones((grad.shape[1],)) # Assuming 2D matrix
	for i in range(grad.shape[1]):	
		beta[i] = \
		residual[:,i].T @ (residual[:,i]-kwargs['residual'][:,i]) / \
		(residual[:,i].T @ residual[:,i])
		ndirection[:,i] = residual[:,i] + beta[i]*kwargs['direction'][:,i]
	return ndirection


def bisection_linmin(atoms, init_step, init_pos, direction, 
	potentials, init_energy):

	pos = np.array(atoms.positions)
	vects = np.array(atoms.get_cell())
	N = len(atoms.positions)
	step = init_step

	while(True):
		# Check if energy is dropping
		pos_temp = init_pos + step*direction
		energy = potentials['Coulomb'].calc(
					pos_array=pos_temp, 
					vects_array=vects, N_=N)['Electrostatic'] + \
				potentials['Buckingham'].calc(
					pos_array=pos_temp, 
					vects_array=vects, N_=N)

		print("Energy difference: ",energy-init_energy)

		if energy>init_energy:
			step = step/2
		else:
			break

		print("New step: ",step)

	# Initialise step
	step_temp = step
	last_energy = energy

	print("Initial step: ",init_step)

	# Try to find better step size
	while(True):
		# Increase step size
		step_temp = step_temp+(step/3)

		print("New step: ",step_temp)

		# Update PES point
		pos_temp = init_pos + step_temp*direction
		energy = potentials['Coulomb'].calc(
					pos_array=pos_temp, 
					vects_array=vects, N_=N)['Electrostatic'] + \
				potentials['Buckingham'].calc(
					pos_array=pos_temp, 
					vects_array=vects, N_=N)

		print("Energy difference: ",energy-init_energy)

		# If new energy is higher abort
		if energy>last_energy:
			return (step, energy)
		elif energy<last_energy:
			step = step_temp
		else:
			raise ValueError("Energy did not change.")

class Descent:
	def __init__(self, ftol=0.00000001, gtol=0.00000001, iterno=1000):
		self.iterno = iterno
		self.ftol = ftol
		self.gtol = gtol

	def iter_step(self, atoms, potentials, step, 
		last_energy, ri=None, pi=None,
		step_func=bisection_linmin, direction_func=GD):
		"""Updating iteration step. The energy is calculated on
		new ion positions that are not inflicted to the Atoms
		object.

		Returns:
		- Direction 	: The direction vector used for the updating step
		- Step 			: The value of step size unchanged
		- Positions 	: The updated ions' positions
		- Gradient 		: The gradient vector w.r.t. the positions before change
		- Energy 		: The value of energy w.r.t. the updated positions
		
		"""
		pos = np.array(atoms.positions)
		vects = np.array(atoms.get_cell())
		N = len(atoms.positions)

		# Gradient
		grad_coul = np.array(
			potentials['Coulomb'].calc_drv(
			pos_array=pos, vects_array=vects, N_=N))
		grad_buck = np.array(
			potentials['Buckingham'].calc_drv(
			pos_array=pos, vects_array=vects, N_=N))
		grad = grad_coul+grad_buck

		# Direction
		pi_1 = direction_func(grad,
			residual=ri, direction=pi)

		# Calculate step size
		(step, energy) = step_func(atoms, step, np.copy(atoms.positions), 
			pi_1, potentials, last_energy)

		# Calculate new point on energy surface
		pos_temp = atoms.positions + step*pi_1

		return {'Direction':pi_1, 'Gradient':grad, 'Step':step, 
		'Positions':pos_temp, 'Energy':energy }

	def first_iteration(self, init_energy, atoms, potentials, 
		step=0.1, step_func=bisection_linmin, direction_func=GD):

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

		p= direction_func(grad,
			residual=-grad, direction=-grad)
		
		# Calculate step size
		(step, energy) = step_func(atoms, step, np.copy(atoms.positions), 
			p, potentials, init_energy)

		# Calculate new point on energy surface
		atoms.positions = atoms.positions + step*p

		gnorm = np.linalg.norm(grad)
		iteration = {'Direction':-grad/gnorm, 'Gradient':grad, 'Step':step, 
		'Positions':atoms.positions, 'Energy':energy}

		prettyprint(iteration)
		input()

		return iteration

	def repeat(self, iteration, atoms, potentials, step=0.1,
		step_func=bisection_linmin, direction_func=GD):

		for i in range(self.iterno):
			last_iteration = iteration
			iteration = self.iter_step(atoms, potentials, step=step,
				last_energy=last_iteration['Energy'],
				ri=-last_iteration['Gradient'], pi=last_iteration['Direction'],
				step_func=step_func, direction_func=direction_func)
			atoms.positions = iteration['Positions']
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