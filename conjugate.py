import sys
import argparse
import fileinput
from ase import *
from ase.io import read
from ase.visualize import view
from ase.calculators.gulp import GULP
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG

''' Get input from file '''
parser = argparse.ArgumentParser(
    description='Define input.')
parser.add_argument(
    'ifilename', metavar='--input', type=str,
    help='Name of input file')
args = parser.parse_args()

print("Input file:" + str(args.ifilename))

''' Set names '''
str = args.ifilename.split("/")
test_no = str[-1].split(".")[0]  # get number of .cif input file
# traj = "CG"+test_no+".traj"
label = "structure"+test_no  # set structure label

''' Read atoms from file '''
atoms = read(args.ifilename)

''' Set calculator. unfix keyword is important; 
by default GULP will set the first derivative 
of the first atom in the list to zero,
since this atom is not usually allowed to move 
during relaxations in GULP '''
calc = GULP(keywords='opti bfgs conp full nosymm',
            library='buck.lib')
atoms.set_calculator(calc)


''' Calculate energy '''
E = atoms.get_potential_energy()
print(E)


''' Perturb atoms '''
# atoms.rattle(stdev=0.1)
# view(atoms)

''' Optimisation with CG internal '''
# calc.set(keywords='conp opti')
# opt = SciPyFminCG(atoms, logfile='-', trajectory=traj,
#                   callback_always=False, alpha=70.0, master=None,
#                   force_consistent=None)
# opt.run()
# E = atoms.get_potential_energy()
# print(E)
