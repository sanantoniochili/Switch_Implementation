import sys
import argparse
import fileinput
from ase import *
from ase.io import read
from ase.visualize import view
from ase.calculators.gulp import GULP

''' Read from file '''
parser = argparse.ArgumentParser(
    description='Define input.')
parser.add_argument(
    'ifilename', metavar='--input', type=str,
    help='Name of input file')
args = parser.parse_args()

print("Reading from file: " + str(args.ifilename))

''' Set calculator and read in previous configurations '''
calc = GULP(keywords='opti bfgs conp full nosymm',
            library='buck.lib')
calc.read(args.ifilename)

''' Calculate energy '''
# E = atoms.get_potential_energy()
# print(E)
