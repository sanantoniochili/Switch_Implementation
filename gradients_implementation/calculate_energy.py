import sys
import shutil
import argparse
import fileinput
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
							   AutoMinorLocator)

from cmath import pi
from cmath import exp
import cmath
import math
from PIL import Image 

from ase import *
from ase.visualize import view
from ase.io import read as aread
from ase.io import write as awrite
from ase.calculators.gulp import GULP
from ase.calculators.lammpslib import LAMMPSlib
from ase.visualize.plot import plot_atoms

from potential import *
from forces import *
from descent import *
from finite_differences import *

DATAPATH = "../../Data/"

charge_dict = {
	'O' : -2.,
	'Sr':  2.,
	'Ti':  4.,
	'Cl': -1.,
	'Na':  1.,
	'S' : -2.,
	'Zn':  2.
}

atom_types = {
	'O':  1,
	'Sr': 2,
	'Ti': 3,
	# 'Na': 4,
	# 'Cl': 5,
	# 'S' : 6,
	# 'Zn': 7
}


def print_template(dict):
	print("\n")
	print("==========================COULOMB==========================".center(columns))
	print("--------------------------CUSTOM IMPLEMENTATION---------------------------------".center(columns))
	print("Real:\t\t"+str(dict['Real']))
	print("Self:\t\t"+str(dict['Self']))
	print("Recip:\t\t"+str(dict['Reciprocal']))
	print("Electrostatic:\t"+str(dict['Electrostatic']))
	print("----------------------------------LAMMPS----------------------------------------".center(columns))
	print("Electrostatic:\t"+str(dict['Elect_LAMMPS']))
	print("---------------------------------MADELUNG---------------------------------------".center(columns))
	if dict['E_madelung'] == None:
		print("No Madelung constant for this structure.")
	else:
		print("Electrostatic:\t"+str(dict['E_madelung']))
	print("--------------------------------------------------------------------------------".center(columns))
	print("\n")
	print("==========================BUCKINGHAM==========================".center(columns))
	print("--------------------------CUSTOM IMPLEMENTATION---------------------------------".center(columns))
	print("Interatomic:\t"+str(dict['Interatomic']))
	print("----------------------------------LAMMPS----------------------------------------".center(columns))
	print("Interatomic:\t"+str(dict['Inter_LAMMPS']))
	print("--------------------------------------------------------------------------------".center(columns))
	print("\n")
	print("==========================TOTAL==========================".center(columns))
	print("--------------------------CUSTOM IMPLEMENTATION---------------------------------".center(columns))
	print("Total lattice:\t"+str(dict['Electrostatic'] + dict['Interatomic']))
	print("-------------------------------------GULP---------------------------------------".center(columns))
	print("Total lattice:\t"+str(dict['Total_GULP']))
	print("--------------------------------------------------------------------------------".center(columns))


def calculate_energies(atoms,accuracy,alpha,real_cut,recip_cut):
	""" Calculate energy using different calculators
	for the Buckingham Coulomb potential. The Coulomb part
	is calculated using the traditional Ewald summation. All
	the parts are printed separately where possible.

	"""
	Er           = 0
	Es 			 = 0
	Erc 		 = 0
	elect_LAMMPS = 0
	Einter 		 = 0
	inter_LAMMPS = 0
	total_GULP 	 = 0
	Emade = None
	LOG = "src/test.log"
	N = len(atoms.positions)
	
	# avoid truncating too many terms
	assert((-np.log(accuracy)/N**(1/6)) >= 1)

	########################### COULOMB ################################

	libfile = DATAPATH+"Libraries/madelung.lib"
	Cpot = Coulomb()
	Cpot.set_parameters(alpha=alpha, real_cut_off=real_cut,
						recip_cut_off=recip_cut, 
						chemical_symbols=atoms.get_chemical_symbols(),
						charge_dict=charge_dict,
						filename=libfile)

	coulomb_energies = Cpot.calc(atoms)

	# # Change atom_style parameter to allow charges
	# LAMMPSlib.default_parameters['lammps_header'] = ['units metal', 'atom_style charge', \
	# 													'atom_modify map array sort 0 0']
	# charge_cmds = [
	# 				"set type 1 charge 2.0",  # Sr
	# 				"set type 2 charge -2.0", # O
	# 				"set type 3 charge 4.0",  # Ti
	# 				# "set type 4 charge 1.0",  # Na
	# 				# "set type 5 charge -1.0", # Cl
	# 				# "set type 6 charge -2.0", # S
	# 				# "set type 7 charge  2.0"  # Zn
	# 				]
	# param_cmds  = [
	# 				"pair_style coul/long "+str(real_cut),
	# 				"pair_coeff * *",
	# 				"kspace_style ewald "+str(accuracy)]
	# cmds = charge_cmds + param_cmds
	# lammps = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, log_file=LOG)
	# atoms.set_initial_charges(Cpot.charges)
	# atoms.set_calculator(lammps)

	# elect_LAMMPS = atoms.get_potential_energy()
	# Emade = Cpot.calc_madelung()

	######################## BUCKINGHAM ################################

	libfile = DATAPATH+"Libraries/buck.lib"
	Bpot = Buckingham()
	Bpot.set_parameters(libfile, atoms.get_chemical_symbols())

	Einter = Bpot.calc(atoms)

	# cmds = [
	# 		"pair_style buck 10",
	# 		"pair_coeff 1 2 1952.39 0.33685 19.22", # Sr-O
	# 		"pair_coeff 1 3 4590.7279 0.261 0.0",   # Ti-O
	# 		"pair_coeff 1 1 1388.77 0.36262 175",   # O-O
	# 		"pair_coeff 2 2 0.0 1.0 0.0",
	# 		"pair_coeff 3 3 0.0 1.0 0.0",
	# 		"pair_coeff 2 3 0.0 1.0 0.0"]
	# lammps = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, log_file=LOG)
	# atoms.set_calculator(lammps)
	# inter_LAMMPS = atoms.get_potential_energy()

	# ############################ TOTAL #################################

	# from ase.calculators.gulp import GULP
	# g_keywords = 'conp full nosymm unfix gradient'
	# calc = GULP(keywords=g_keywords, \
	# 				library='buck.lib')
	# atoms.set_calculator(calc)
	# total_GULP = atoms.get_potential_energy()

	# dict = { **coulomb_energies,
	# 		'Elect_LAMMPS': elect_LAMMPS, 'E_madelung': Emade, 'Interatomic': Einter,
	# 		'Inter_LAMMPS': inter_LAMMPS, 'Total_GULP': total_GULP}
	# print_template(dict)

	return {'Coulomb': Cpot, 'Buckingham': Bpot, \
			'energy': coulomb_energies['Electrostatic']+Einter}


def calculate_forces(atoms, potentials):
	"""Calculation of forces using the analytical
	derivative of the energy.

	"""
	pos = atoms.positions
	vects = atoms.get_cell()
	N = len(pos)

	dcoul = DCoulomb(potentials['Coulomb'])
	grad_real = dcoul.calc_real(pos, vects, N)
	grad_recip = dcoul.calc_recip(pos,vects, N)
	grad_coul = grad_real + grad_recip

	# dbuck = DBuckingham(potentials['Buckingham'])
	# grad_buck = dbuck.calc(pos, vects, N)
	# grad = grad_coul+grad_buck
	grad = grad_coul
	gnorm = np.linalg.norm(grad)
	return {'grad': grad, 'gnorm': gnorm}


def get_input(filename=None):
	"""Choose input from file or manual.

	"""
	folder = None
	structure = None
	if filename:
		atoms = aread(filename)
		structure = "structure"+args.i.rstrip('.cif').split('/')[-1]
		folder = "random"
		print("Using file as input.")
	else:
		atoms = Atoms("SrTiO3",

				  cell=[[4.00, 0.00, 0.00],
						[0.00, 4.00, 0.00],
						[0.00, 0.00, 4.00]],

				  positions=[[0, 0, 0],
							 [2, 2, 2],
							 [0, 2, 2],
							 [2, 0, 2],
							 [2, 2, 0]],
				  pbc=True)
		# atoms = atoms.repeat((3, 1, 1))
		print("Using custom Atoms object as input.")
	return (folder,structure,atoms)


if __name__ == "__main__":
	columns = shutil.get_terminal_size().columns
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'-i', metavar='--input', type=str,
		help='.cif file to read')
	args = parser.parse_args()
	(folder, structure, atoms) = get_input(args.i)


	######################### ENERGY ################################
	''' parameters '''
	vects = np.array(atoms.get_cell())
	volume = abs(np.linalg.det(vects))
	N = len(atoms.get_positions())
	accuracy = 0.00001  # Demanded accuracy of terms 
	alpha = N**(1/6) * pi**(1/2) / volume**(1/3)
	real_cut = (-np.log(accuracy))**(1/2)/alpha
	recip_cut = 2*alpha*(-np.log(accuracy))**(1/2)

	''' calculation '''
	potentials = calculate_energies(\
		atoms,accuracy,alpha,real_cut,recip_cut)
	print(atoms.get_chemical_symbols())
	print("---Initial atom positions")
	print(atoms.positions)
	# view(atoms)

	######################### FIN DIFFS ##############################
	# np.set_printoptions(precision=10,suppress=True)
	diffs = finite_diff_grad(atoms=atoms, ions=range(N), coords=range(3), \
							 initial_energy=potentials['energy'],\
							 displacement=0.01, Coulomb=potentials['Coulomb'],\
							 Buckingham=potentials['Buckingham'], N=N, vects=vects)
	print("---Finite differences:")
	print(diffs)
	diffs_ = finite_diff_grad(atoms=atoms, ions=range(N), coords=range(3), \
							 initial_energy=potentials['energy'],\
							 displacement=-0.01, Coulomb=potentials['Coulomb'],\
							 Buckingham=potentials['Buckingham'], N=N, vects=vects)
	print("----> f(x+h)")
	print(diffs)
	print("----> f(x-h)")
	print(diffs_)
	print("----> f(x+h)-f(x-h)")
	print(abs(diffs_)-abs(diffs))

	######################### FORCES #################################
	forces = calculate_forces(atoms, potentials)['grad']
	print("---Analytical gradient")
	print(forces)

	atoms.positions[0] += [0.01,0,0]
	print("---New atom positions")
	print(atoms.positions)

	forces_ = calculate_forces(atoms, potentials)['grad']
	print("---New analytical gradient")
	print(forces_)

	print("---Difference of initial gradient vs moved ion gradient")
	print(forces-forces_)
	
	######################### ANGLE #################################
	# print(np.round(forces['grad'],decimals=10))
	# fforces = flatten(forces)
	# fdiffs = flatten(diffs)
	# acos = np.dot(fforces,fdiffs) / (np.linalg.norm(fforces)*np.linalg.norm(fdiffs))
	# print("---Diffs and analytical gradient angle:")
	# print(acos)

	# print(forces-diffs)

	######################### RELAXATION #############################
	# GDescent = Descent()
	# GDescent.repeat(atoms, potentials, 0.1)


# https://github.com/SINGROUP/Pysic/blob/master/fortran/Geometry.f90
# https://github.com/vlgusev/IPCSP/blob/master/tools/matrix_generator.py?

def flatten(positions):
	vector = []
	for sublist in positions:
		for element in sublist:
			vector.append(element)
	return vector