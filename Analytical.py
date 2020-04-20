import sys
import argparse
import fileinput
import numpy as np
from scipy.special import erfc

from cmath import pi
from cmath import exp
import cmath
import math

from ase import *
from ase.visualize import view
from ase.io import read as aread
from ase.geometry import Cell

from BCpot import Potential

DATAPATH = "../../Data/"

charge_dict = {
	'O' : -2.,
	'Sr':  2.,
	'Ti':  4.,
	'Cl': -1.,
	'Na':  1.
}

if __name__=="__main__":
	atoms  = aread(DATAPATH+"material/NaCl.cif")
	print(atoms.get_positions())
	vects  = np.array(atoms.get_cell())
	volume = abs(np.linalg.det(vects))
	alpha  = 2/(volume**(1.0/3))

	################ COULOMB ################
	cutoff = 4
	pot = Potential()
	pot.set_structure(charge_dict, atoms)
	shifts = pot.get_shifts( cutoff,pot.vects )

	esum = 0
	for ioni in range(0, pot.N): ### wrong if N not 2
		for ionj in range(ioni, pot.N): 
			dist = pot.pos[ionj,] - pot.pos[ioni,]
			if ioni != ionj:
				esum += ( 1/np.linalg.norm(dist) )
			for shift in shifts:
				esum += ( -1/(np.linalg.norm(dist+shift)) )

	print(esum)

