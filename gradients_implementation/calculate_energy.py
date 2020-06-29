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
import output_functions as outf

DATAPATH = "../../../Data/"

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


def calculate_energies(atoms):
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

	######################## INITIALISATION ############################

	vects = np.array(atoms.get_cell())
	volume = abs(np.linalg.det(vects))
	N = len(atoms.get_positions())
	accuracy = 0.00001  # Demanded accuracy of terms 
	alpha = N**(1/6) * pi**(1/2) / volume**(1/3)
	real_cut = (-np.log(accuracy))**(1/2)/alpha
	recip_cut = 2*alpha*(-np.log(accuracy))**(1/2)
	LOG = "src/test.log"

	if -np.log(accuracy)/N**(1/6) < 1:
		raise ValueError("Values of constants truncate too many terms.")

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

	dict = { **coulomb_energies,
			'Elect_LAMMPS': elect_LAMMPS, 'E_madelung': Emade, 'Interatomic': Einter,
			'Inter_LAMMPS': inter_LAMMPS, 'Total_GULP': total_GULP}
	print_template(dict)

	return {'Coulomb': Cpot, 'Buckingham': Bpot, \
			'energy': coulomb_energies['Electrostatic']+Einter}


def calculate_forces(atoms, potentials):
	"""Calculation of forces using the analytical
	derivative of the energy.

	"""
	dcoul = DCoulomb(potentials['Coulomb'])
	grad_real = dcoul.calc_real(atoms)
	grad_recip = dcoul.calc_recip(atoms)
	grad_coul = grad_real + grad_recip

	dbuck = DBuckingham(potentials['Buckingham'])
	grad_buck = dbuck.calc(atoms)
	grad = grad_coul+grad_buck

	gnorm = np.linalg.norm(grad)
	return {'grad': grad, 'gnorm': gnorm}


def finite_diff_grad(atoms, initial_energy, displacement):
	"""Defining local slope using finite differences 
	(like the derivative limit). 

	A displacement is added to the coordinate of each ion per iter.
	The potential <new_energy> is evaluated using the new 
	position and then we have the partial derivative's
	approach as (<new_energy> - <initial_energy>) / <displacement>
	for each one of x iterations (trials).

	"""
	dims = 3
	grad = np.zeros((len(atoms.positions),3))
	h = np.zeros(3)
	for ioni in range(len(atoms.positions)):
		for coord in range(dims):
			h[:dims] = 0
			h[coord] = displacement
			# print("Ion: {} moved {}".format(ioni, h))
			# careful not to exceed unit cell
			vects_coords = [vect[coord] for vect in atoms.get_cell()]
			# try not to move ioni out of unit cell
			if atoms.positions[ioni][coord]+h[coord] <= max(vects_coords):
				# add perturbation to one ion coordinate
				atoms.positions[ioni] += h
				# calculate (f(x+h)-f(x))/h
				grad[ioni][coord] = ( calculate_energies(atoms)['energy']-initial_energy )/h[coord]
				# restore initial coordinates
				atoms.positions[ioni] -= h
			else:
				# add perturbation to one ion coordinate
				atoms.positions[ioni] -= h
				# calculate (f(x+h)-f(x))/h
				grad[ioni][coord] = ( calculate_energies(atoms)['energy']-initial_energy )/h[coord]
				# restore initial coordinates
				atoms.positions[ioni] += h
				print("Ion: {} moved the other way".format(ioni))
	return grad


def get_input(filename):
	"""Choose input from file or manual.

	"""
	if filename:
		atoms = aread(filename)
		structure = "structure"+args.i.rstrip('.cif').split('/')[-1]
		folder = "random"
	else:
		atoms = Atoms("NaCl",

	              cell=[[2.00, 0.00, 0.00],
	                    [0.00, 2.00, 0.00],
	                    [0.00, 0.00, 2.00]],

	              positions=[[0, 0, 0],
	                         [1, 1, 1],
	                         ],
	              pbc=True)
		atoms = atoms.repeat((3, 1, 1))
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


	displacement = 0.01
	''' ENERGY '''
	potentials = calculate_energies(atoms)

	''' FORCES '''
	# forces = calculate_forces(atoms, potentials)
	# diffs = finite_diff_grad(atoms, potentials['energy'], displacement)
	# acos = np.dot(forces['grad'],diffs)/(forces['gnorm']*np.linalg.norm(diffs))

	# print("Diffs vs Forces angle cos: {}".format(acos))
	# outf.print_forces_vs_diffs(folder, structure, forces['grad'], diffs, displacement)

	''' RELAXATION '''
	# GDescent = Descent()
	# GDescent.repeat(atoms, potentials)

	# for x in range(len(atoms.get_positions())):
	# 	acos = np.dot(forces['grad'][x],diffs[x])
	# 	print(acos)

# https://github.com/SINGROUP/Pysic/blob/master/fortran/Geometry.f90
# https://github.com/vlgusev/IPCSP/blob/master/tools/matrix_generator.py?
