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


DATAPATH = "../../Data/"


charge_dict = {
	'O': -2.,
	'Sr': 2.,
	'Ti': 4.}

class Potential:
	def __init__(self):
		pass

	def set_structure(self, charge_dict, atoms):
		'''
		 atoms:       ASE object
		 vects:       Unit cell vectors
		 positions:   ASE positions of atoms
		 charge_dict: Dictionary of chemical symbols with charge value
		'''
		self.vects   = atoms.get_cell()
		self.pos     = atoms.get_positions()
		self.volume  = abs(np.linalg.det(self.vects))
		self.charges = [charge_dict[x] for x in atoms.get_chemical_symbols()]
		self.N       = len(self.pos)
		self.atoms   = atoms

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


class Coulomb(Potential):
	'''
	 Calculations for the Coulomb energy contribution.
	 Ewald summation method used for long range.
	'''
	def __init__(self, alpha, real_cut_off=4, recip_cut_off=4):
		self.real_cut_off  = real_cut_off
		self.recip_cut_off = recip_cut_off
		self.alpha         = alpha

	def get_reciprocal_vects(self):
		'''
		 Calculate reciprocal vectors
		'''
		rvects = np.zeros((3,3))
		for i in np.arange(3):
			rvects[i,] = 2*pi*np.cross( self.vects[(1+i)%3,], \
									self.vects[(2+i)%3] ) / self.volume
		return rvects

	def get_charges_mult(self, index1, index2):
		'''
		 Find respective charges and
		 return their product
		'''
		return ( self.charges[index1]*self.charges[index2] )

	def calc_self(self, esum=[]):
		'''
		 Calculate self interaction term
		'''
		if esum == []:
			esum = np.zeros((Cpot.N, Cpot.N))
		for i in range(0, self.N):
			esum[i, i] -= ( self.get_charges_mult(i, i) * ( alpha / math.sqrt(pi) ))
		return esum

	def calc_real(self, esum=[]):
		'''
		 Calculate short range
		'''
		if esum == []:
			esum = np.zeros((Cpot.N, Cpot.N))
		shifts = self.get_shifts( self.real_cut_off,self.vects )
		for ioni in range(0, self.N): 
			for ionj in range(ioni, self.N): 
				if ioni != ionj: # skip in case it's the same atom in original unit cell
					rij = np.linalg.norm(self.pos[ioni,] - self.pos[ionj,])
					esum[ioni, ionj] += ( self.get_charges_mult(ioni,ionj) * \
												 math.erfc( self.alpha*rij )/(2*rij) )
				# take care of the rest lattice (+ Ln)
				for shift in shifts:
					rij = np.linalg.norm(self.pos[ioni,] + shift - self.pos[ionj,])
					esum[ioni, ionj] += ( self.get_charges_mult(ioni,ionj) * \
												math.erfc( self.alpha*rij )/(2*rij) )
		return esum

	def calc_recip(self, recip_vects, esum=[]):
		'''
		 Calculate long range
		'''
		if esum == []:
			esum        = np.zeros((Cpot.N, Cpot.N))
		shifts = self.get_shifts( self.recip_cut_off,recip_vects )
		for ioni in range(0, self.N): 
			for ionj in range(ioni, self.N): 

				dist = self.pos[ionj,] - self.pos[ioni,] 
				for k in shifts:
					po = -np.dot(k,k)/(4*alpha**2)
					numerator = 4 * (pi**2) * (math.exp(po)) * math.cos(np.dot(k, dist))
					denominator = np.dot(k,k) * 2 * pi * self.volume
					esum[ioni, ionj] += (( self.get_charges_mult(ioni,ionj) ) * \
																(numerator/denominator))
		return esum

	def calc_complete(self, esum):
		'''
		 Complete lower triangular matrix of monopole to monopole
		'''
		esum *= 14.399645351950543
		for ioni in range(0, self.N):
			for ionj in range(0, ioni):
				esum[ioni, ionj] = esum[ionj, ioni]
		return esum

class Buckingham(Potential):
	'''
	 Calculations for the Buckingham energy contribution.
	'''
	def __init__(self, filename, real_cut_off=2):
		'''
		 Set atom_i-atom_j parameters as read from library file:
		 - par: [A(eV), rho(Angstrom), C(eVAngstrom^6)]
		 - lo : min radius (Angstrom)
		 - hi : max radius (Angstrom)
		'''
		self.buck = {}
		self.real_cut_off  = real_cut_off

		try:
			with open(filename,"r") as fin:
				for line in fin:
					line = line.split()
					if ( len(line)<4 ):
						continue
					pair = (min(line[0],line[2]), max(line[0],line[2]))
					self.buck[pair] = {}
					self.buck[pair]['par'] = list(map(float, line[4:7]))
					self.buck[pair]['lo'] = float(line[7])
					self.buck[pair]['hi'] = float(line[-1])
		except IOError:
			print("No library file found.")

	def calc_real(self, esum=[]):
		'''
		 Interatomic potential for original unit cell
		'''
		if esum == []:
			esum        = np.zeros((Cpot.N, Cpot.N))
		chemical_symbols = self.atoms.get_chemical_symbols()
		for ioni in range(self.N):
			for ionj in range(ioni, self.N):
				# Find the pair we are examining
				pair = (min(chemical_symbols[ioni], chemical_symbols[ionj]), \
								max(chemical_symbols[ioni], chemical_symbols[ionj]))
				if (pair in self.buck):
				# Pair of ions is listed in parameters file
					A    = self.buck[pair]['par'][0]
					rho  = self.buck[pair]['par'][1]
					C    = self.buck[pair]['par'][2]

					dist = np.linalg.norm(self.pos[ioni] - self.pos[ionj])
					# Check if distance of ions allows interaction 					
					if (dist <= self.buck[pair]['hi']) & (ioni != ionj):
						esum[ioni, ionj] += A*math.exp(-1.0*dist/rho) - C/dist**6
					
					# Check interactions with neighbouring cells
					shifts = self.get_shifts( int(self.buck[pair]['hi']),self.vects )
					for shift in shifts:
						dist = np.linalg.norm(self.pos[ioni] + shift - self.pos[ionj])
						# Check if distance of ions allows interaction 					
						if (dist <= self.buck[pair]['hi']):
							esum[ioni, ionj] += A*math.exp(-1.0*dist/rho) - C/dist**6
		return esum

	def calc_complete(self, esum):
		'''
		 Complete lower triangular matrix of monopole to monopole
		'''
		for ioni in range(0, self.N):
			for ionj in range(0, ioni+1):
				esum[ioni, ionj] = esum[ionj, ioni]
		return esum

if __name__=="__main__":
	atoms  = aread(DATAPATH+"RandomStart_Sr3Ti3O9/3.cif")
	vects  = np.array(atoms.get_cell())

	volume = abs(np.linalg.det(vects))
	alpha  = 2/(volume**(1.0/3))

	Cpot        = Coulomb(alpha,4,4)
	Cpot.set_structure(charge_dict, atoms)
	rvects      = Cpot.get_reciprocal_vects()

	Er  = Cpot.calc_real()
	Es  = Cpot.calc_self()
	Erc = Cpot.calc_recip(rvects)

	Eupper = Er + Es + Erc
	Etotal = Cpot.calc_complete(Eupper)

	print("--------------------------------------------------------------------------------")

	# Convert to eV per Angstrom
	print("Real:\t\t"+str(sum(sum(Cpot.calc_complete(Er)))))
	print("Self:\t\t"+str(sum(sum(Cpot.calc_complete(Es)))))
	print("Recip:\t\t"+str(sum(sum(Cpot.calc_complete(Erc)))))
	print("Total:\t\t"+str(sum(sum(Etotal))))

	print("--------------------------------------------------------------------------------")

	filename    = DATAPATH+"Libraries/buck.lib"
	Bpot 		= Buckingham(filename)
	Bpot.set_structure(charge_dict, atoms)

	Einter = Bpot.calc_real()
	print("Interatomic:\t"+str(sum(sum(Einter))))

	print("--------------------------------------------------------------------------------")

	print("Total lattice:\t"+str(sum(sum(Etotal+Einter))))

	print("--------------------------------------------------------------------------------")

	# from matrix_generator import BuckinghamTwoIons
	# for ioni in range(Bpot.N):
	# 	for ionj in range(ioni+1, Bpot.N):
	# 		# Find the pair we are examining
	# 		pair = (min(chemical_symbols[ioni], chemical_symbols[ionj]), \
	# 						max(chemical_symbols[ioni], chemical_symbols[ionj]))
	# 		if (pair in self.buck):
	# 		# Pair of ions is listed in parameters file
	# 			A    = self.buck[pair]['par'][0]
	# 			rho  = self.buck[pair]['par'][1]
	# 			C    = self.buck[pair]['par'][2]

	# 			dist = np.linalg.norm(self.pos[ioni] - self.pos[ionj])
	# 			# Check if distance of ions allows interaction 					
	# 			if (dist <= self.buck[pair]['hi']):
	# 				esum[ioni, ionj] += A*math.exp(-1.0*dist/rho) - C/dist**6
				
	# 			# Check interactions with neighbouring cells
	# 			shifts = self.get_shifts( int(self.buck[pair]['hi']),self.vects )
	# 			for shift in shifts:
	# 				dist = np.linalg.norm(self.pos[ioni] + shift - self.pos[ionj])
	# 				# Check if distance of ions allows interaction 					
	# 				if (dist <= self.buck[pair]['hi']):
	# 					esum[ioni, ionj] += A*math.exp(-1.0*dist/rho) - C/dist**6

	# print("Inter:\t"+str(Einter))
	# print(sum(sum(Einter)))

# https://github.com/SINGROUP/Pysic/blob/master/fortran/Geometry.f90
# https://github.com/vlgusev/IPCSP/blob/master/tools/matrix_generator.py?
