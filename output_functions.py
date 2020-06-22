import os,sys
import shutil
import argparse
import fileinput
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

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
	
