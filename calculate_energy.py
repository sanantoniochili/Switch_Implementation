import sys
import argparse
import fileinput
import numpy as np
import pandas as pd

from cmath import pi
from cmath import exp
import cmath
import math

from ase import *
from ase.visualize import view
from ase.io import read as aread
from ase.calculators.gulp import GULP
from ase.calculators.lammpslib import LAMMPSlib

from potential import *
from forces import *

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
	'O' : 1,
	'Sr': 2,
	'Ti': 3,
	'Na': 4,
	'Cl': 5,
	'S' : 6,
	'Zn': 7
}

import shutil
if __name__=="__main__":
	columns = shutil.get_terminal_size().columns
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'filename', metavar='--input', type=str,
		help='.cif file to read')
	args = parser.parse_args()
	atoms = aread(args.filename)


	############################################################################################
	###################################### INITIALISATION ######################################
	############################################################################################

	vects    = np.array(atoms.get_cell())
	volume   = abs(np.linalg.det(vects))
	N        = len(atoms.get_positions())
	accuracy = 0.00001 # Demanded accuracy of terms (tolerance of parameters?)
	alpha    = N**(1/6) * pi**(1/2) / volume**(1/3)
	LOG      = "src/test.log"

	#########################################################################################
	######################################## COULOMB ########################################
	#########################################################################################
	
	print("\n")
	print("==========================COULOMB==========================".center(columns))

	libfile    = DATAPATH+"Libraries/madelung.lib"
	real_cut  = (-np.log(accuracy))**(1/2)/alpha
	recip_cut = 2*alpha*(-np.log(accuracy))**(1/2)
	Cpot      = Coulomb(charge_dict, atoms)
	Cpot.set_parameters(alpha,real_cut_off=real_cut,recip_cut_off=recip_cut, filename=libfile)

	print("--------------------------CUSTOM IMPLEMENTATION---------------------------------".center(columns))
	rvects    = Cpot.get_reciprocal_vects()
	Er  = Cpot.calc_real()
	Es  = Cpot.calc_self()
	Erc = Cpot.calc_recip(rvects)
	Eupper = Er + Es + Erc
	Etotal = Cpot.calc_complete(Eupper)

	print("Real:\t\t"+str(sum(sum(Cpot.calc_complete(Er)))))
	print("Self:\t\t"+str(sum(sum(Cpot.calc_complete(Es)))))
	print("Recip:\t\t"+str(sum(sum(Cpot.calc_complete(Erc)))))
	print("Electrostatic:\t"+str(sum(sum(Etotal))))
	print("----------------------------------LAMMPS----------------------------------------".center(columns))

	# Change atom_style parameter to allow charges
	LAMMPSlib.default_parameters['lammps_header'] = ['units metal', 'atom_style charge', \
														'atom_modify map array sort 0 0']
	charge_cmds = [
					"set type 1 charge 2.0",  # Sr
					"set type 2 charge -2.0", # O
					"set type 3 charge 4.0",  # Ti
					"set type 4 charge 1.0",  # Na
					"set type 5 charge -1.0", # Cl
					"set type 6 charge -2.0", # S
					"set type 7 charge  2.0"  # Zn
					]
	param_cmds  = [
					"pair_style coul/long "+str(real_cut),
					"pair_coeff * *",
					"kspace_style ewald "+str(accuracy)]
	cmds = charge_cmds + param_cmds
	lammps = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, log_file=LOG)
	atoms.set_initial_charges(Cpot.charges)
	atoms.set_calculator(lammps)

	print("Electrostatic:\t"+str(atoms.get_potential_energy()))
	print("---------------------------------MADELUNG---------------------------------------".center(columns))
	Emade = Cpot.calc_madelung()
	if Emade==None:
		print("No Madelung constant for this structure.")
	else:
		print("Electrostatic:\t"+str(Emade))
	print("--------------------------------------------------------------------------------".center(columns))

	############################################################################################
	######################################## BUCKINGHAM ########################################
	############################################################################################

	# print("\n")
	# print("==========================BUCKINGHAM==========================".center(columns))

	# print("--------------------------CUSTOM IMPLEMENTATION---------------------------------".center(columns))
	libfile    = DATAPATH+"Libraries/buck.lib"
	Bpot 	   = Buckingham(charge_dict, atoms)
	Bpot.set_parameters(libfile)
	Einter = Bpot.calc()
	# print("Interatomic:\t"+str(Einter))
	# print("----------------------------------LAMMPS----------------------------------------".center(columns))
	
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
	
	# print("Interatomic:\t"+str(atoms.get_potential_energy()))
	# print("--------------------------------------------------------------------------------".center(columns))

	############################################################################################
	########################################### TOTAL ##########################################
	############################################################################################

	print("\n")
	print("==========================TOTAL==========================".center(columns))
	print("--------------------------CUSTOM IMPLEMENTATION---------------------------------".center(columns))
	print("Total lattice:\t"+str(sum(sum(Etotal)) + Einter))

	print("-------------------------------------GULP---------------------------------------".center(columns))
	from ase.calculators.gulp import GULP
	calc = GULP(keywords='conp full nosymm unfix', \
					library='buck.lib')
	atoms.set_calculator(calc)
	print("Total lattice:\t"+str(atoms.get_potential_energy()))

	print("--------------------------------------------------------------------------------".center(columns))

	############################################################################################
	########################################## OUTPUT ##########################################
	############################################################################################

	# ''' SAVE TO CSV '''
	# ''' NAME '''
	# count = args.filename.split('.')[-2].split('/')[-1]
	# structure = 'structure'+count
	# ''' FOLDER '''
	# folder = 'random'
	# ''' DATAFRAME '''
	# df = pd.DataFrame.from_dict({'structure': [structure], 'folder': [folder], 'interatomic_custom': [Einter],\
	# 							 'interatomic_lammps': [atoms.get_potential_energy()], 'electrostatic_custom': [sum(sum(Etotal))] }, orient='columns')
	# fileout = "src/custom_all.csv"
	# try:
	# 	df_in = pd.read_csv(fileout) 
	# 	if(df_in.empty):
	# 		df.to_csv(fileout, mode='w', header=True)
	# 	else:
	# 		df.to_csv(fileout, mode='a', header=False)
	# except:
	# 	df.to_csv(fileout, mode='w', header=True)

	# dcoul = DCoulomb(Cpot)
	# grad_real = dcoul.calc_real()
	# grad_recip = dcoul.calc_recip()	
	# grad_coul = grad_real + grad_recip

	# dbuck = DBuckingham(Bpot)
	# grad_buck = dbuck.calc()
	# grad = grad_coul+grad_buck

	# gnorm = np.linalg.norm(grad) 

	# print(gnorm)
	


# https://github.com/SINGROUP/Pysic/blob/master/fortran/Geometry.f90
# https://github.com/vlgusev/IPCSP/blob/master/tools/matrix_generator.py?
