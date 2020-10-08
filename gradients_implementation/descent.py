import argparse
import numpy as np
import pandas as pd

from ase.io import read as aread
from cysrc.potential import *

class Descent:
	def __init__(self):
		pass

	def calculate_direction(self, atoms, potentials):
		pos = atoms.positions
		vects = atoms.get_cell()
		N = len(pos)

		dcoul = DCoulomb(potentials['Coulomb'])
		grad_coul = dcoul.calc_real(pos, vects, N)
		grad_coul += dcoul.calc_recip(pos , vects, N)

		dbuck = DBuckingham(potentials['Buckingham'])
		grad_buck = dbuck.calc(pos, vects, N)
		grad = grad_coul+grad_buck

		return -grad

	def repeat(self, atoms, potentials):
		x_energy = potentials['Coulomb'].calc(atoms)['Electrostatic'] + \
							potentials['Buckingham'].calc(atoms)
		"""Iteration step.

		"""
		count = 0
		for x in range(100):
			# Direction
			p = self.calculate_direction(atoms, potentials)
			print("Iter: {} \tEnergy: {} \tDirection: {} \tStep: {}".format(\
				count,x_energy,np.linalg.norm(p),step),flush=True)

			# Step
			step = 0.5

			# Calculate new point on energy surface
			pos_temp = np.copy(atoms.positions + step*p)

			# Calculate new energy 
			vects = atoms.get_cell()
			N = len(atoms.positions)
			energy = \
			potentials['Coulomb'].calc(\
				None, positions=pos_temp, vects=vects, N=N)['Electrostatic'] + \
							potentials['Buckingham'].calc_real(pos_temp, vects, N)

			# If new energy is lower, keep it
			if energy<x_energy:
				atoms.positions = pos_temp
				x_energy = energy
				count += 1
				continue

			step = step/4
			count += 1
		print("Iter: {} \tEnergy: {} \tDirection: {} \tStep: {}".format(\
				count,energy,np.linalg.norm(p),step),flush=True)
			


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