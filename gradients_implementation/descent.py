import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from cysrc.potential import *

def GD(atoms, potentials, grad, **kwargs):
	"""Gradient Descent updating scheme. Returns
	the direction vector.

	"""
	return -grad/np.linalg.norm(grad)

def CG(atoms, potentials, grad, **kwargs):
	"""Conjugate Gradient updating scheme. We are
	using Polak-Ribière updating for beta factor. 
	Returns the direction vector.

	"""

	residual = -grad
	beta = np.dot(residual,(residual-kwargs['residual']))/np.dot(residual,residual)
	return residual+beta*kwargs['direction']


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
	p = direction_func(atoms, potentials, grad,
		residual=ri, direction=pi)

	# Calculate new point on energy surface
	pos_temp = np.copy(atoms.positions + step*p)

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

	p= direction_func(atoms, potentials, grad,
		residual=-grad, direction=-grad)

	# Calculate new point on energy surface
	atoms.positions = atoms.positions + step*p

	# Calculate new energy 
	energy = potentials['Coulomb'].calc(
			pos_array=pos, 
			vects_array=vects, N_=N)['Electrostatic'] + \
		potentials['Buckingham'].calc(
			pos_array=pos, 
			vects_array=vects, N_=N)

	iteration = {'Direction':-grad, 'Gradient':grad, 'Step':step, 
	'Positions':pos, 'Energy':energy}

	# Rest iterations
	for i in range(10):
		last_iteration = iteration
		iteration = iter_step(atoms, potentials, step,
			-last_iteration['Gradient'], last_iteration['Direction'])
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