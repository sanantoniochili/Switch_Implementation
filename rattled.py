from ase import *
import sys
import argparse
from ase.io import read, write
from ase.visualize import view
from ase.calculators.gulp import GULP

''' Get stdev from user '''
parser = argparse.ArgumentParser(
    description='Define input.')
parser.add_argument(
    'stdev', metavar='--input', type=float,
    help='Name of input file')
args = parser.parse_args()

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


''' Perturb atoms '''
atoms.rattle(stdev=args.stdev)

''' Wriet to .cif files '''
filename = "input/rattled/dev"+str(args.stdev)+".cif"
write(filename, atoms, format='cif')