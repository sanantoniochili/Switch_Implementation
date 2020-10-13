import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from cysrc.potential import *

def GD(potentials, grad, **kwargs):
	"""Gradient Descent updating scheme. Returns
	the direction vector.

	"""
	return -grad/np.linalg.norm(grad)

def CG(potentials, grad, **kwargs):
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


def iter_step(atoms, potentials, step, 
	ri=None, pi=None, direction_func=GD):
	"""Updating iteration step. The energy is calculated on
	new ion positions that are not inflicted to the Atoms
	object.

	Returns:
	- Direction 	: The direction vector used for the updating step
	- Step 			: The value passed to the function unchanged
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
	pi_1 = direction_func(potentials, grad,
		residual=ri, direction=pi)

	# Calculate new point on energy surface
	pos_temp = np.copy(atoms.positions + step*pi_1)

	# Calculate new energy 
	energy = potentials['Coulomb'].calc(
			pos_array=pos_temp, 
			vects_array=vects, N_=N)['Electrostatic'] + \
		potentials['Buckingham'].calc(
			pos_array=pos_temp, 
			vects_array=vects, N_=N)
	return {'Direction':p, 'Gradient':grad, 'Step':step, 
	'Positions':pos_temp, 'Energy':energy }


def repeat(atoms, potentials, direction_func=GD):

	step = 0.1

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

	p= direction_func(potentials, grad,
		residual=-grad, direction=-grad)
	return
	# Calculate new point on energy surface
	atoms.positions = atoms.positions + step*p

	# Calculate new energy 
	energy = potentials['Coulomb'].calc(
			pos_array=pos, 
			vects_array=vects, N_=N)['Electrostatic'] + \
		potentials['Buckingham'].calc(
			pos_array=pos, 
			vects_array=vects, N_=N)

	gnorm = np.linalg.norm(grad)
	iteration = {'Direction':-grad/gnorm, 'Gradient':grad, 'Step':step, 
	'Positions':atoms.positions, 'Energy':energy}

	print(iteration)

	# Rest iterations
	for i in range(10):
		last_iteration = iteration
		iteration = iter_step(atoms, potentials, step,
			ri=-last_iteration['Gradient'], pi=last_iteration['Direction'],
			direction_func=direction_func)
		atoms.positions = iteration['Positions']
		print(iteration)

			


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'-i', metavar='--input', type=str,
		help='.cif file to read')
	args = parser.parse_args()
	atoms = aread(args.i)