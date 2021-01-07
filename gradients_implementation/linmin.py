import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from cysrc.potential import *
from ase.geometry import get_distances
from math import log
import sys


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
	max_step = -log(sys.float_info.max)/100

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
		# print("Initial energy: ",init_iter['Energy'],
		# 	" New energy: ",energy, "Step: ",step, "\nDirection: ",end="")
		# for x in direction:
		# 	print(x, end="")
		# print("\n")

		# Vectorise matrices
		vdirection = np.reshape(direction, 
			(direction.shape[0]*direction.shape[1],1))
		vgrad = np.reshape(grad, 
			(grad.shape[0]*grad.shape[1],1))

		if (energy-init_iter['Energy']) <= step*c1*(vdirection.T @ vgrad):
			break

	return (step, energy)