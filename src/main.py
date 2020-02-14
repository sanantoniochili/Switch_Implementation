# import torch
import sys
import argparse
import fileinput

from cmath import pi
from cmath import exp
import cmath

from ase import *
from ase.io import read as aread
from ase.geometry import Cell

library = {
	'O'  : -2.,
	'Sr' : 2.,
	'Ti' : 4.}
		

Atoms = aread("/users/phd/tonyts/Desktop/Data/RandomStart_Sr3Ti3O9/1.cif")
charges =  [library[x] for x in Atoms.get_chemical_symbols()]
Atoms.set_initial_charges(charges)

''' Buckingham-Coulomb '''

rs = Atoms.get_positions()
rCell = 2*pi*Atoms.get_reciprocal_cell()  # 2pi not included in ase
vol = Atoms.get_volume()

'''
 - phi = potential field 
 - e0  = vacuum permitivity
'''

e0 = 8.854 * 10**(-12)
sigma = 0.1
El = 0 # long distance energy

# for all nonzero recipr vects
for ks in rCell:
	for k in ks:
		if not ( len(rs)==len(charges) ): # something wrong with atom number
			raise ValueError
		Sk = 0
		for atom in range(0,len(rs)):
			# Sk += charges[atom]*exp( -----k*rs[atom] -- dot product ----- *1j)
		El += (exp(- sigma**2 * k**2)/2) * (abs(Sk)**2) 


print(El)