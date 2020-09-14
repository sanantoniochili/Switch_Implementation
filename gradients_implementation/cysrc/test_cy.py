from ase import *
from cmath import pi

import sys
import numpy as np
sys.path.append('/home/sanantoniochili/Desktop/PhD/Scripts/Switch_Implementation/gradients_implementation/pysrc')
from example import Potential, Coulomb, Buckingham
from potential import Coulomb as PCoulomb
from potential import Buckingham as BP

DATAPATH = "/home/sanantoniochili/Desktop/PhD/Data/"

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

charge_dict = {
	'O' : -2.,
	'Sr':  2.,
	'Ti':  4.,
	'Cl': -1.,
	'Na':  1.,
	'S' : -2.,
	'Zn':  2.
}


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


if __name__=="__main__":
	import shutil, argparse
	columns = shutil.get_terminal_size().columns
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'-i', metavar='--input', type=str,
		help='.cif file to read')
	args = parser.parse_args()
	(folder, structure, atoms) = get_input(args.i)

	vects = np.array(atoms.get_cell())
	volume = abs(np.linalg.det(vects))
	N = len(atoms.get_positions())
	accuracy = 0.00001  # Demanded accuracy of terms 
	alpha = N**(1/6) * pi**(1/2) / volume**(1/3)
	real_cut = (-np.log(accuracy))**(1/2)/alpha
	recip_cut = 2*alpha*(-np.log(accuracy))**(1/2)

	p = Coulomb()
	pp = PCoulomb()

	p.set_parameters(alpha=alpha, real_cut_off=real_cut,
							recip_cut_off=recip_cut, 
							chemical_symbols=np.array(atoms.get_chemical_symbols()),
							charge_dict=charge_dict, N=N,
							filename=None)
	pp.set_parameters(alpha=alpha, real_cut_off=real_cut,
							recip_cut_off=recip_cut, 
							chemical_symbols=atoms.get_chemical_symbols(),
							charge_dict=charge_dict,
							filename=None)
	# pos = atoms.positions
	# pos = np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1],[1,1,1,1]], dtype=np.double)
	# vects=np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1]], dtype=np.double)


	libfile = DATAPATH+"Libraries/buck.lib"
	Bpot = Buckingham()
	Bpot.set_parameters(libfile, np.array(atoms.get_chemical_symbols()))

	PBpot = BP()
	PBpot.set_parameters(libfile, np.array(atoms.get_chemical_symbols()))

	import timeit
	pt = timeit.timeit('''print(Bpot.calc(atoms)+p.calc(atoms)['Electrostatic'])''', globals=globals(), number=1)
	ppt = timeit.timeit('''print(PBpot.calc(atoms)+pp.calc(atoms)['Electrostatic'])''', globals=globals(), number=1)

	print('Cython: %f' % pt)
	print('Python: %f' % ppt)