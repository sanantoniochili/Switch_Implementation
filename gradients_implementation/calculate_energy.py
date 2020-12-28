import sys
import shutil
import argparse
import fileinput
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

from cysrc.potential import *
from finite_differences import finite_diff_grad

import timeit

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
	import shutil
	columns = shutil.get_terminal_size().columns
	print("\n")

	# Coulomb
	lines = columns/2 - len("COULOMB")
	for c in range(lines):
		print("=", end="")
	print("COULOMB", end="")
	for c in range(lines):
		print("=", end="")
	print()

	lines = columns/2 - len("CUSTOM IMPLEMENTATION")
	for c in range(lines):
		print("-", end="")
	print("CUSTOM IMPLEMENTATION", end="")
	for c in range(lines):
		print("-", end="")
	print()

	print("Real:\t\t"+str(dict['Real']))
	print("Self:\t\t"+str(dict['Self']))
	print("Recip:\t\t"+str(dict['Reciprocal']))
	print("Electrostatic:\t"+str(dict['Electrostatic']))

	lines = columns/2 - len("LAMMPS")
	for c in range(lines):
		print("-", end="")
	print("LAMMPS", end="")
	for c in range(lines):
		print("-", end="")
	print()

	print("Electrostatic:\t"+str(dict['Elect_LAMMPS']))

	lines = columns/2 - len("MADELUNG")
	for c in range(lines):
		print("-", end="")
	print("MADELUNG", end="")
	for c in range(lines):
		print("-", end="")
	print()

	if dict['E_madelung'] == None:
		print("No Madelung constant for this structure.")
	else:
		print("Electrostatic:\t"+str(dict['E_madelung']))

	for c in range(columns):
		print("-", end="")
	print()
	print("\n")
	
	# Buckingham
	lines = columns/2 - len("BUCKINGHAM")
	for c in range(lines):
		print("=", end="")
	print("BUCKINGHAM", end="")
	for c in range(lines):
		print("=", end="")
	print()

	lines = columns/2 - len("CUSTOM IMPLEMENTATION")
	for c in range(lines):
		print("-", end="")
	print("CUSTOM IMPLEMENTATION", end="")
	for c in range(lines):
		print("-", end="")
	print()

	print("Interatomic:\t"+str(dict['Interatomic']))

	lines = columns/2 - len("LAMMPS")
	for c in range(lines):
		print("-", end="")
	print("LAMMPS", end="")
	for c in range(lines):
		print("-", end="")
	print()

	print("Interatomic:\t"+str(dict['Inter_LAMMPS']))
	
	for c in range(columns):
		print("-", end="")
	print()
	print("\n")

	# Total
	lines = columns/2 - len("TOTAL")
	for c in range(lines):
		print("=", end="")
	print("TOTAL", end="")
	for c in range(lines):
		print("=", end="")
	print()

	lines = columns/2 - len("CUSTOM IMPLEMENTATION")
	for c in range(lines):
		print("-", end="")
	print("CUSTOM IMPLEMENTATION", end="")
	for c in range(lines):
		print("-", end="")
	print()

	print("Total lattice:\t"+str(dict['Electrostatic'] + dict['Interatomic']))
	
	lines = columns/2 - len("GULP")
	for c in range(lines):
		print("-", end="")
	print("GULP", end="")
	for c in range(lines):
		print("-", end="")
	print()
	
	print("Total lattice:\t"+str(dict['Total_GULP']))
	
	for c in range(columns):
		print("-", end="")
	print()


def lammps_energy(atoms):
	LAMMPSlib.default_parameters['lammps_header'] = ['units metal', 'atom_style charge', \
														'atom_modify map array sort 0 0']
	charge_cmds = [
					"set type 1 charge 2.0",  # Sr
					"set type 2 charge -2.0", # O
					"set type 3 charge 4.0",  # Ti
					# "set type 4 charge 1.0",  # Na
					# "set type 5 charge -1.0", # Cl
					# "set type 6 charge -2.0", # S
					# "set type 7 charge  2.0"  # Zn
					]
	param_cmds  = [
					"pair_style coul/long "+str(real_cut),
					"pair_coeff * *",
					"kspace_style ewald "+str(accuracy)]
	cmds = charge_cmds + param_cmds
	lammps = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, log_file=LOG)
	atoms.set_initial_charges(Cpot.charges)
	atoms.set_calculator(lammps)
	elect_LAMMPS = atoms.get_potential_energy()

	cmds = [
			"pair_style buck 10",
			"pair_coeff 1 2 1952.39 0.33685 19.22", # Sr-O
			"pair_coeff 1 3 4590.7279 0.261 0.0",   # Ti-O
			"pair_coeff 1 1 1388.77 0.36262 175",   # O-O
			"pair_coeff 2 2 0.0 1.0 0.0",
			"pair_coeff 3 3 0.0 1.0 0.0",
			"pair_coeff 2 3 0.0 1.0 0.0"]
	lammps = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, log_file=LOG)
	atoms.set_calculator(lammps)
	inter_LAMMPS = atoms.get_potential_energy()

	return (elect_LAMMPS,inter_LAMMPS)


def gulp_energy(atoms):
	from ase.calculators.gulp import GULP
	g_keywords = 'conp full nosymm unfix gradient'
	calc = GULP(keywords=g_keywords, \
					library='buck.lib')
	atoms.set_calculator(calc)
	total_GULP = atoms.get_potential_energy()
	return total_GULP


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
							 # [1.5, .5, 2], # put Os too close
							 [2, 2, 0]],
				  pbc=True)
		# atoms = atoms.repeat((3, 1, 1))
		print("Using custom Atoms object as input.")
		view(atoms)
	return (folder,structure,atoms)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'-i', metavar='--input', type=str,
		help='.cif file to read')
	parser.add_argument(
		'-r', '--relax', action='store_true',
		help='Perform structural relaxation')
	parser.add_argument(
		'-rl', '--lagrangian', action='store_true',
		help='Use Lagrangian relaxation')
	args = parser.parse_args()
	(folder, structure, atoms) = get_input(args.i)

	####################### INITIALISE ##############################

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

	######################### ENERGY ################################
	
	''' parameters '''
	vects = np.array(atoms.get_cell())
	volume = abs(np.linalg.det(vects))
	N = len(atoms.get_positions())
	accuracy = 0.00001  # Demanded accuracy of terms 
	alpha = N**(1/6) * pi**(1/2) / volume**(1/3)
	real_cut = (-np.log(accuracy))**(1/2)/alpha
	recip_cut = 2*alpha*(-np.log(accuracy))**(1/2)
	chemical_symbols=np.array(atoms.get_chemical_symbols())

	# avoid truncating too many terms
	assert((-np.log(accuracy)/N**(1/6)) >= 1)

	libfile = DATAPATH+"Libraries/madelung.lib"
	Cpot = Coulomb()
	Cpot.set_parameters(alpha=alpha, real_cut_off=real_cut,
						recip_cut_off=recip_cut, 
						chemical_symbols=chemical_symbols,
						N=N, charge_dict=charge_dict,
						filename=libfile)
	coulomb_energies = Cpot.calc(atoms)
	
	libfile = DATAPATH+"Libraries/buck.lib"
	libfile2 = DATAPATH+"Libraries/radii.lib"
	Bpot = Buckingham()
	Bpot.set_parameters(filename=libfile, 
		chemical_symbols=chemical_symbols, radius_lib=libfile2)
	Einter = Bpot.calc(atoms)

	# dict = { **coulomb_energies,
	# 		'Elect_LAMMPS': elect_LAMMPS, 'E_madelung': Emade, 'Interatomic': Einter,
	# 		'Inter_LAMMPS': inter_LAMMPS, 'Total_GULP': total_GULP}
	# print_template(dict)

	# ######################## TIMING #################################

	# import cProfile
	# import re
	# import pstats
	# from pstats import SortKey
	# cProfile.run(\
	# 	'calculate_energies(atoms,accuracy,alpha,real_cut,recip_cut)',\
	# 	 'profiling_energy.log')
	# p = pstats.Stats('profiling_energy.log')
	# print("\nStats for potentials:")
	# print("Cumulative stats:")
	# p.sort_stats(SortKey.CUMULATIVE).print_stats(11)
	# print("Total time stats:")
	# p.sort_stats(SortKey.TIME).print_stats(10)

	# cProfile.run('calculate_forces(atoms, potentials)', 'profiling_dervs.log')
	# p = pstats.Stats('profiling_dervs.log')
	# print("Stats for derivatives:")
	# print("Cumulative stats:")	
	# p.sort_stats(SortKey.CUMULATIVE).print_stats(11)
	# print("Total time stats:")
	# p.sort_stats(SortKey.TIME).print_stats(10)

	######################### RELAXATION #############################
	from descent import *

	desc = Descent()
	initial_energy = coulomb_energies['Electrostatic']+Einter
	potentials = {'Coulomb':Cpot, 'Buckingham':Bpot}

	# print(atoms.positions)
	# print(chemical_symbols)
	# Cpot.calc_real_drv2(np.array(atoms.positions), vects, N)

	# print("ALL DISTS:")
	# print(atoms.get_distances())
	Bpot.catastrophe_check(atoms.positions, 0.5)

	if args.relax:
		if args.lagrangian:
			libfile = DATAPATH+"Libraries/radii.lib"
			LCpot = Lagrangian()
			potentials =  {**potentials, 'Lagrangian':LCpot}

			vlambdas = np.random.rand(int(N*(N-1)/2))
			assert(np.all(vlambdas>0))
			lambdas = np.zeros((N,N))
			count = vlambdas.shape[0]

			# Populate lambda matrix with a value for
			# each pairwise distance

			lambdas[0][4], lambdas[4][0] = .25,.25
			lambdas[2][3], lambdas[3][2] = .1,.1

			print("Printing lambda vector..")
			for ioni in range(N):
				print()
				for ionj in range(ioni+1, N):
					# lambdas[ioni][ionj] = vlambdas[count-1]
					# lambdas[ionj][ioni] = vlambdas[count-1]
					print("Ions {},{} : {}".format(\
						ioni, ionj, lambdas[ioni][ionj]), end="  ")
					count -= 1
			print()

			input()

			print("Using lambda multipliers vector: {}".format(lambdas))
			LCpot.set_parameters(lambdas, libfile, chemical_symbols)
			initial_energy += LCpot.calc_constrain(
				pos=atoms.positions, 
				N=N)
		
		iters = desc.repeat(
			init_energy=initial_energy,
			atoms=atoms, 
			potentials=potentials, 
			direction_func=CG,
			step_func=bisection_linmin)

	