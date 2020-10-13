import os,sys
import shutil
import argparse
import fileinput
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.visualize import view
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
# import pysrc.potential as pt

def stepsize_energy_change(foldr, filename, atoms):
	"""Move ion with different stepsize to check 
	changes in energy.

	"""
	steps = []
	energies = []

	step = 0.01
	ioni = 0
	coord = 0
	for x in range(11):
		atoms.positions[ioni] += [step,0,0]
		energies += [calculate_energy(atoms)['energy']]
		steps += [x*step]

	fileout = "src/displacement_echange/"+folder+"_"+structure+"_echange"

	df = pd.DataFrame({'steps':steps, 'energies':energies})
	df.to_csv(fileout+".csv")
	ax = df.plot('steps', 'energies', figsize=(10,10))
	ax.ticklabel_format(useOffset=False)

	plt.xlabel('steps')
	plt.ylabel('energy')
	plt.savefig(fileout+".png")
	plt.show()


def finite_diff_grad(atoms, ions, coords, initial_energy, displacement, **kwargs):
	"""Defining local slope using finite differences 
	(like the derivative limit). 

	A displacement is added to the coordinate of each ion per iter.
	The potential <new_energy> is evaluated using the new 
	position and then we have the partial derivative's
	approach as (<new_energy> - <initial_energy>) / <displacement>
	for each one of x iterations (trials).

	"""
	dims = len(coords)
	grad = np.zeros((len(atoms.positions),3))
	h = np.zeros(3)

	for ioni in ions:
		for coord in coords:
			h[:dims] = 0
			h[coord] = displacement
			# try not to move ioni out of unit cell 
			# in the affected direction
			vects_coords = [vect[coord] for vect in atoms.get_cell()]
			assert(abs(atoms.positions[ioni][coord]+h[coord]) \
												<= max(vects_coords))
			# add perturbation to one ion coordinate
			init_pos = atoms.positions[ioni][coord]
			positions_cp = atoms.positions.copy()
			positions_cp[ioni] += h
			new_pos = positions_cp[ioni][coord]

			# print("\n-- Ion {}{} moved {} in dimension {}".format(\
			# 	atoms.get_chemical_symbols()[ioni],ioni,displacement,coord))
			# print("-- {} --> {}".format(\
			# 	init_pos,new_pos))

			pos = np.array(positions_cp)
			vects = np.array(atoms.get_cell())
			N = len(atoms.positions)
			# calculate (f(x+h)-f(x))/h
			grad[ioni][coord] = ( kwargs['Buckingham'].calc(None, pos, vects, N)+\
				kwargs['Coulomb'].calc(None, pos, vects, N)['Electrostatic']\
				-initial_energy )/h[coord]

	return grad


def print_forces_vs_diffs(folder, struct, forces, diffs, displacement):
	"""Print dataframe with finite differences moving every ion to 
	every direction and analytical derivative.

	"""
	df = pd.DataFrame(diffs)
	df.columns = ['diffs_' + str(col) for col in df.columns]
	df['displacement'] = displacement
	df = df.set_index('ion_'+df.index.astype(str))

	dff = pd.DataFrame(forces)
	dff.columns = ['forces_' + str(col) for col in dff.columns]
	dff = dff.set_index('ion_'+dff.index.astype(str))
	df = df.join(dff)

	arrays = [[folder], [struct], df.index.tolist()]
	index = pd.MultiIndex.from_product(arrays, names=('folder', 'structure', 'ion'))
	df.set_index(index, inplace=True)
	print(df)

	fileout = "src/finite_differences/fin_diffs_all_coords_all_atoms_3_structs.csv"
	if not os.path.isfile(fileout):
		try:
			with open(fileout, 'w') as f:
				df.to_csv(f, header=True)
		finally:
			f.close()
	else:
		try:
			with open(fileout, 'a') as f:
				df.to_csv(f, header=False)
		finally:
			f.close()
	
def forces_diffs_angle(folder, struct, forces, diffs):
	diffs = flatten(diffs)
	forces = flatten(forces)

	df = pd.DataFrame(diffs)
	df = df.T.join(pd.DataFrame(forces).T, lsuffix='diffs', rsuffix='forces')
	df['cos'] = np.dot(diffs,forces)/(np.linalg.norm(diffs)*np.linalg.norm(forces))
	print("Diffs vs Forces angle cos: {}".format(np.dot(diffs,forces)/(np.linalg.norm(diffs)*np.linalg.norm(forces))))

	arrays = [[folder], [struct]]
	index = pd.MultiIndex.from_product(arrays, names=('folder', 'structure'))
	df.set_index(index, inplace=True)

	# fileout = "finite_differences/diffs_forces_angles.csv"
	# if not os.path.isfile(fileout):
	# 	try:
	# 		with open(fileout, 'w') as f:
	# 			df.to_csv(f, header=True)
	# 	finally:
	# 		f.close()
	# else:
	# 	try:
	# 		with open(fileout, 'a') as f:
	# 			df.to_csv(f, header=False)
	# 	finally:
	# 		f.close()

def flatten(positions):
	vector = []
	for sublist in positions:
		for element in sublist:
			vector.append(element)
	return vector