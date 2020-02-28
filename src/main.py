# import torch
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
from ase.io import read as aread
from ase.geometry import Cell


DATAPATH = "/users/phd/tonyts/Desktop/Data/"


charge_dict = {
	'O': -2.,
	'Sr': 2.,
	'Ti': 4.}

class Supercell:
	'''
	 Create supercell
	'''
	def __init__(self, size, atoms):
		'''
		 Unit cell possible positions are size^3 combinations
		 Every unit cell has N atoms
		 Every atoms has a position of 3 coordinates
		'''
		self.atoms     = atoms
		self.size      = size
		self.vects     = [[size,0,0],[0,size,0],[0,0,size]]@atoms.get_cell()
		self.positions = np.zeros(((2*size+1)**3, len(atoms.get_positions()), 3))

	def populate(self):
		'''
		 Create atoms supercell positions based on unit cell
		'''
		i = 0
		size_all_directions = 2*self.size+1
		# For all possible cell positions (combinations of integers)
		for spot in np.ndindex(size_all_directions, size_all_directions, \
													size_all_directions):
			# Cell position
			[shift_x,shift_y,shift_z] = spot-np.array([self.size,self.size,self.size])
			# Multiply cell coords with respactive atom coords
			self.positions[i] = [ [pos[0]*shift_x, pos[1]*shift_y, pos[2]*shift_z]  \
											for pos in self.atoms.get_scaled_positions()]
			i += 1

class Coulomb:
	'''
	 Calculations for the Coulomb energy contribution.
	 Ewald summation method used for long range.
	'''
	def __init__(self, alpha, real_cut_off=4, recip_cut_off=4):
		self.real_cut_off  = real_cut_off
		self.recip_cut_off = recip_cut_off
		self.alpha         = alpha

	def set_structure(self, charge_dict, scell, volume):
		'''
		 atoms:       ASE object
		 vects:       Supercell vectors
		 pos:         ASE positions of atoms
		 charge_dict: Dictionary of chemical symbols with charge value
		'''
		self.vects   = scell.vects
		self.pos     = scell.positions
		self.volume  = volume
		# self.charges = [charge_dict[x] for x in atoms.get_chemical_symbols()]
		self.N       = len(scell.positions)

		print("--------------------------------------------------------------")
		print("Vectors used: ")
		print(self.vects)
		print("Vectors of unit cell:")
		print(atoms.get_cell()[:])
		print("\n")
		print("Atoms positions:")
		print([list(p) for p in self.pos])
		print("\n")
		print("Volume:")
		print(volume)
		print("\n")
		print("Atoms charges:")
		print(self.charges)
		print("--------------------------------------------------------------")

	def get_reciprocal_vects(self):
		'''
		 Calculate reciprocal vectors
		'''
		rvects = np.zeros((3,3))
		for i in np.arange(3):
			rvects[i,] = 2*pi*np.cross( self.vects[(1+i)%3,], \
									self.vects[(2+i)%3] ) / self.volume
		return rvects

	def get_shifts(self, cut_off, vects):
		''' 
		 Get all possible lattice positions:   
		 (2cut_off+1)^3-{case of (cut_off,cut_off,cut_off)}
		 combinations in R^3  
		'''
		shifts = np.zeros(((2*cut_off+1)**3 -1, 3))
		tmp = np.array([cut_off,cut_off,cut_off])

		i = 0
		for shift in np.ndindex(2*cut_off+1, 2*cut_off+1, 2*cut_off+1):
			if shift != (cut_off,cut_off,cut_off):
				shifts[i,] = shift
				shifts[i,] = shifts[i,] - tmp
				i = i+1
		shifts = shifts@vects
		return shifts

	def get_charges_mult(self, index1, index2):
		'''
		 Find respective charges and
		 return their product
		'''
		return ( self.charges[index1]*self.charges[index2] )

	def calc_self(self): # CHECK CHARGES AND POSITIONS
		'''
		 Calculate self interaction term
		'''
		esum        = np.zeros((Cpot.N, Cpot.N))
		for i in range(0, self.N):
			esum[i, i] += ( self.get_charges_mult(i, i) )
			esum[i, i] *= ( alpha / math.sqrt(pi) )
		return esum

# 	def calc_real(self):
# 		'''
# 		 Calculate short range
# 		'''
# 		esum        = np.zeros((Cpot.N, Cpot.N))
# 		shifts = self.get_shifts( self.real_cut_off,self.vects )
# 		for atomi in range(0, self.N): 
# 			for atomj in range(atomi, self.N): 

# 				if atomi != atomj: # skip in case it's the same atom
# 					rij = np.linalg.norm(self.pos[atomi,] - self.pos[atomj,])
# 					esum[atomi, atomj] += math.erfc( self.alpha*rij )/(2*rij)

# 					# take care of the rest lattice (+ Ln)
# 					for shift in shifts: # !!!!!!! MAYBE NEEDS TO BE OUTSIDE if BLOCK !!!!!!!
# 						rij = np.linalg.norm(self.pos[atomi,] + shift - self.pos[atomj,])
# 						esum[atomi, atomj] += math.erfc( self.alpha*rij )/(2*rij)

# 					esum[atomi, atomj] *= ( self.get_charges_mult(atomi,atomj) )
# 		return esum

# 	def calc_recip(self, recip_vects):
# 		'''
# 		 Calculate long range
# 		'''
# 		esum        = np.zeros((Cpot.N, Cpot.N))
# 		shifts = self.get_shifts( self.recip_cut_off,recip_vects )
# 		for atomi in range(0, self.N): 
# 			for atomj in range(atomi, self.N): 

# 				dist = self.pos[atomj,] - self.pos[atomi,] # order of atomi, atomj matters

# 				for k in shifts:
# 					# po = -np.dot(k,k)/(4*alpha**2)
# 					# numerator = 4 * (pi**2) * (math.exp(po)) * math.cos(np.dot(k, dist))
# 					# denominator = np.dot(k,k) * 2 * pi * self.volume
# 					# esum[atomi, atomj] += (( self.get_charges_mult(atomi,atomj) ) * \
# 					# 											(numerator/denominator))
# 					term = (4*math.pi**2)/np.dot(k,k)
# 					term = term*math.exp(-np.dot(k,k)/(4*alpha**2))
# 					term = term*math.cos(np.dot(k, dist))
# 					esum[atomi, atomj] += term/(2*math.pi*self.volume)

# 				esum[atomi, atomj] *= ( self.get_charges_mult(atomi,atomj) )
# 		return esum

# 	def calc_complete(self, esum):
# 		'''
# 		 Complete lower triangular matrix of monopole to monopole
# 		'''
# 		for atomi in range(0, self.N):
# 			for atomj in range(0, atomi):
# 				esum[atomi, atomj] = esum[atomj, atomi]
# 		return esum


# class Buckingham:
# 	'''
# 	 Calculations for the Buckingham energy contribution.
# 	'''
# 	def __init__(self, filename):
# 		'''
# 		 Set atom_i-atom_j parameters as read from library file:
# 		 - par: [A(eV), rho(Angstrom), C(eVAngstrom^6)]
# 		 - lo : min radius (Angstrom)
# 		 - hi : max radius (Angstrom)
# 		'''
# 		self.buck = {}

# 		try:
# 			with open(filename,"r") as fin:
# 				for line in fin:
# 					line = line.split()
# 					if ( len(line)<4 ):
# 						continue
# 					pair = (min(line[0],line[2]), max(line[0],line[2]))
# 					self.buck[pair] = {}
# 					self.buck[pair]['par'] = list(map(float, line[4:7]))
# 					self.buck[pair]['lo'] = float(line[7])
# 					self.buck[pair]['hi'] = float(line[-1])
# 		except IOError:
# 			print("No library file found.")


if __name__=="__main__":
	atoms = aread("/users/phd/tonyts/Desktop/Data/RandomStart_Sr3Ti3O9/1.cif")
	
	size         = 2
	scell        = Supercell(size, atoms)
	scell.populate()

	cell_volume = abs(np.linalg.det(scell.vects))
	alpha       = 2/(cell_volume**(1.0/3))

	Cpot        = Coulomb(alpha,4,4)
	# Cpot.set_structure(scell.positions, cell_volume, charge_dict)
	# rvects      = Cpot.get_reciprocal_vects()

	# Er = Cpot.calc_real()
	# # print(Er)	
	# Es = Cpot.calc_self()
	# # print(Es)
	# Er = Cpot.calc_recip(rvects)
	# # print(Er)

	# Eupper = Er + Es + Er
	# Etotal = Cpot.calc_complete(Eupper)

	# # Convert to eV per Angstrom
	# Etotal *= 14.399645351950543
	# print(Etotal.min())
	# print(Etotal.max())

	# filename    = DATAPATH+"Libraries/buck.lib"
	# Bpot 		= Buckingham(filename)

# https://github.com/SINGROUP/Pysic/blob/master/fortran/Geometry.f90
# https://github.com/vlgusev/IPCSP/blob/master/tools/matrix_generator.py?

