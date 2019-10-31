import sys
import argparse
import fileinput
from ase import *
from ase.io import read
from ase.visualize import view
from ase.calculators.gulp import GULP

''' Get input from file '''
parser = argparse.ArgumentParser(
    description='Define input.')
parser.add_argument(
    'input', metavar='--input', type=str,
    help='Name of input file')
args = parser.parse_args()

print("Input file:" + str(args.input))

''' Read atoms from file '''
atoms = read(args.input)

''' Set calculator. unfix keyword is important; 
by default GULP will set the first derivative 
of the first atom in the list to zero,
since this atom is not usually allowed to move 
during relaxations in GULP '''
calc = GULP(keywords='gradient conp unfix full nosymm',
            library='buck.lib')
atoms.calc = calc

''' Calculate energy '''
E = atoms.get_potential_energy()
print(E)

''' Optimisation with GULP's internal '''
calc.set(keywords='conp opti')
opt = calc.get_optimizer(atoms)
opt.run()
E = atoms.get_potential_energy()
print(E)
