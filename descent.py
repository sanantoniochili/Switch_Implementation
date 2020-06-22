import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread

from potential import *
from forces import *

class Descent:
	def __init__(self):
		pass

	def calculate_direction(self, potentials):
		dcoul = DCoulomb(potentials['Coulomb'])
		grad_coul = dcoul.calc_real()
		grad_coul += dcoul.calc_recip()

		dbuck = DBuckingham(potentials['Buckingham'])
		grad_buck = dbuck.calc()
		grad = grad_coul+grad_buck

		return -grad

	def repeat(self, potentials):
		# for x in range(10):
		p = self.calculate_direction(potentials)
		# .atoms.positions = self.atoms.positions + p
			


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description='Define input')
	parser.add_argument(
		'-i', metavar='--input', type=str,
		help='.cif file to read')
	args = parser.parse_args()
	atoms = aread(args.i)

	descent = Descent(atoms)
	view(descent.atoms)