import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from cysrc.potential import *
from ase.geometry import get_distances

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


def interpolation_linmin(atoms, direction, potentials, grad, 
	init_iter, c1=0, c2=1):

	vects = np.array(atoms.get_cell())
	N = len(atoms.positions)
	step0 = 1 

	# Check if energy is dropping
	pos_temp = atoms.positions + step0*direction
	energy = potentials['Coulomb'].calc(
				pos_array=pos_temp, 
				vects_array=vects, N_=N)['Electrostatic'] + \
			potentials['Buckingham'].calc(
				pos_array=pos_temp, 
				vects_array=vects, N_=N)
	if "Lagrangian" in potentials:
		energy += potentials['Lagrangian'].calc_constrain(
			pos=pos_temp, N=N)

	# Vectorise matrices
	vdirection = np.reshape(direction, 
		(direction.shape[0]*direction.shape[1],1))
	vgrad = np.reshape(grad, 
		(grad.shape[0]*grad.shape[1],1))

	# Make sure that the new step is large enough to decrease energy
	if (energy-init_iter['Energy']) <= step0*c1*(vdirection.T @ vgrad):
		return (step0, energy)

	# New step size
	step1 = - step0**2 * (vgrad.T @ vdirection) / \
	2 * (energy - init_iter['Energy'] - (vgrad.T @ vdirection)*step0)
	
	# Check if energy is dropping
	pos_temp = atoms.positions + step1*direction
	energy1 = potentials['Coulomb'].calc(
				pos_array=pos_temp, 
				vects_array=vects, N_=N)['Electrostatic'] + \
			potentials['Buckingham'].calc(
				pos_array=pos_temp, 
				vects_array=vects, N_=N)
	if "Lagrangian" in potentials:
		energy1 += potentials['Lagrangian'].calc_constrain(
			pos=pos_temp, N=N)

	# Make sure that the new step is large enough to decrease energy
	if (energy1-init_iter['Energy']) <= step1*c1*(vdirection.T @ vgrad):
		return (step1, energy1)

	energy0 = init_iter['Energy']
	aa_mult = 1/(step0**2 * step1**2 * (step1 - step0))
	phi_grad = vdirection.T @ vgrad

	# Reserve for cubic interpolation
	phi0 = energy
	phi1 = energy1
	stepi = step0
	stepi_1 = step1

	while(True):
	
		# Cubic interpolation parameters
		alpha = aa_mult*(stepi**2 * \
			(phi1-energy0-stepi_1*phi_grad) - \
			stepi_1**2 * (phi0-energy0-stepi*phi_grad))
		beta = aa_mult*(-stepi**3*(phi1-energy0-stepi_1*phi_grad) + \
			stepi_1**3 * (phi0-energy0-stepi*phi_grad))
		print(beta**2-3*alpha*phi_grad)
		nstep = (-beta+(beta**2-3*alpha*phi_grad)**(1/2))/(3*alpha)

		# Check if energy is dropping
		pos_temp = atoms.positions + nstep*direction
		nenergy = potentials['Coulomb'].calc(
					pos_array=pos_temp, 
					vects_array=vects, N_=N)['Electrostatic'] + \
				potentials['Buckingham'].calc(
					pos_array=pos_temp, 
					vects_array=vects, N_=N)
		if "Lagrangian" in potentials:
			nenergy += potentials['Lagrangian'].calc_constrain(
				pos=pos_temp, N=N)

		# Make sure that the new step is large enough to decrease energy
		if (nenergy-init_iter['Energy']) <= nstep*c1*(vdirection.T @ vgrad):
			return (nstep, nenergy)

		# Keep last two calculated values of energy and step size
		phi0 = phi1
		phi1 = nenergy
		stepi = stepi_1
		stepi_1 = nstep


def bisection_linmin(atoms, direction, potentials, grad, 
	init_iter, c1=0, min_step=0.000000001):

	vects = np.array(atoms.get_cell())
	N = len(atoms.positions)
	max_step = 1/min_step

	energy = init_iter['Energy']
	step = 1

	print("In line search:")
	while(step<=max_step):
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

		# Vectorise matrices
		vdirection = np.reshape(direction, 
			(direction.shape[0]*direction.shape[1],1))
		vgrad = np.reshape(grad, 
			(grad.shape[0]*grad.shape[1],1))
		
		# Wolfe conditions
		# If new energy is large then try to decrease step size
		if (energy-init_iter['Energy']) > step*c1*(vdirection.T @ vgrad):
			break

		# Increase step size while acceptable
		step = step*2
		print("New step: ",step)
	
	while(step>=min_step):
		# Decrease step size until accepted
		step = step/2
		print("New step: ",step)

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

		# Vectorise matrices
		vdirection = np.reshape(direction, 
			(direction.shape[0]*direction.shape[1],1))
		vgrad = np.reshape(grad, 
			(grad.shape[0]*grad.shape[1],1))

		if (energy-init_iter['Energy']) <= step*c1*(vdirection.T @ vgrad):
			break

	return (step, energy)


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

	def check_completion(self, init_energy, energy, init_gnorn=None, gnorm=None):
		de = abs(init_energy-energy)
		if de < self.ftol:
			print("Iterations: {} Energy: {} Final Energy difference: {} Tolerance: {}".format(
					self.iters,energy,de,self.ftol))
			return True
		return False

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
			last_iter, self.ftol)

		# Calculate new point on energy surface
		self.iters += 1
		pos_temp = atoms.positions + step*pi_1

		return {'Positions':pos_temp, 'Direction':pi_1, 'Gradient':grad, 
		'Iter':self.iters, 'Step':step, 'Energy':energy}

	def repeat(self, init_energy, atoms, potentials, 
		step_func=bisection_linmin, direction_func=GD):

		pos = np.array(atoms.positions)
		vects = np.array(atoms.get_cell())
		N = len(atoms.positions)

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
			'Direction': p}, self.ftol)

		# Calculate new point on energy surface
		self.iters += 1
		atoms.positions = atoms.positions + step*p
		# Catastrophe check
		if potentials['Buckingham'].catastrophe_check(\
			atoms.positions, self.cattol) is not None:
			print("!! Detected ions too close !!".center(COLUMNS))

		iteration = {'Positions':atoms.positions, 'Direction':p, 
		'Gradient':grad, 'Iter':self.iters, 'Step':step,  'Energy':energy}

		prettyprint(iteration)
		input()

		if self.check_completion(init_energy, iteration['Energy']):
			return iteration

		# Rest of iterations ##########################################
		for i in range(self.iterno):
			last_iteration = iteration
			iteration = self.iter_step(atoms, potentials, 
				last_iter=last_iteration, step_func=step_func, 
				direction_func=direction_func)
			
			# Keep the newly found energy value
			self.emin = iteration['Energy']
			# Change point on PES and try a whole step
			atoms.positions = iteration['Positions']
			# Catastrophe check
			if potentials['Buckingham'].catastrophe_check(\
				atoms.positions, self.cattol) is not None:
				print("!! Detected ions too close !!".center(COLUMNS))

			if self.check_completion(last_iteration['Energy'], iteration['Energy']):
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