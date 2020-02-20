from ase import *
import sys
import argparse
import numpy as np
from ase.io import read, write
from ase.visualize import view
from ase.calculators.gulp import GULP

''' Get stdev from user '''
parser = argparse.ArgumentParser(
    description='Define input.')
parser.add_argument(
    'stdev', metavar='--input', type=float,
    help='Standard deviation')
parser.add_argument(
    'draw', metavar='--draw',
    help='Number of same derivation trial')
args = parser.parse_args()

# print("stdev: "+str(args.stdev))
print("draw: "+str(args.draw))

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
atoms = atoms.repeat((3, 1, 1))

''' Perturb atoms '''
atoms.rattle(stdev=args.stdev, rng=np.random.RandomState(np.random.seed()))

''' Write to .cif files '''
filename = "../../Data/rattled/dev"+str(args.stdev)+"/dev"+str(args.stdev)+"draw"+args.draw+".cif"
write(filename, atoms, format='cif')
