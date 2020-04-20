from scipy.special import erfc
import numpy as np

from cmath import pi
from cmath import exp
import cmath
import math

from ase import *
from ase.visualize import view
from ase.geometry import Cell

class Potential:
	def __init__(self, charge_dict, atoms):
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

	def set_parameters(self):
		pass

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
	def set_parameters(self, alpha, real_cut_off=4, recip_cut_off=4, filename=None):
		self.real_cut_off  = math.ceil(real_cut_off)
		self.recip_cut_off = math.ceil(recip_cut_off)
		self.alpha         = alpha
		self.made_const    = 0

		try:
			with open(filename,"r") as fin:
				for line in fin:
					line = line.split()
					found = True
					for symbol in self.atoms.get_chemical_symbols():
						if symbol not in line: # chemical formula does not match
							found = False
							break
					if found == True: # correct formula
						self.made_const = float(line[-1])						
		except IOError:
			print("No Madelung library file found.")

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
			esum = np.zeros((self.N, self.N))
		for i in range(0, self.N):
			esum[i, i] -= ( self.get_charges_mult(i, i) * ( self.alpha / math.sqrt(pi) ))
		return esum

	def calc_real(self, esum=[]):
		'''
		 Calculate short range
		'''
		if esum == []:
			esum = np.zeros((self.N, self.N))
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
			esum        = np.zeros((self.N, self.N))
		shifts = self.get_shifts( self.recip_cut_off,recip_vects )
		for ioni in range(0, self.N): 
			for ionj in range(ioni, self.N): 

				dist = self.pos[ionj,] - self.pos[ioni,] 
				for k in shifts:
					po = -np.dot(k,k)/(4*self.alpha**2)
					numerator = 4 * (pi**2) * (math.exp(po)) * math.cos(np.dot(k, dist))
					denominator = np.dot(k,k) * 2 * pi * self.volume
					esum[ioni, ionj] += (( self.get_charges_mult(ioni,ionj) ) * \
																(numerator/denominator))
		return esum

	def calc_complete(self, esum):
		'''
		 Complete lower triangular matrix of monopole to monopole
		'''
		esum *= 14.399645351950543 # electrostatic constant
		for ioni in range(0, self.N):
			for ionj in range(0, ioni):
				esum[ioni, ionj] = esum[ionj, ioni]
		return esum

	def calc_madelung(self):
		if not self.made_const:
			return None
		esum = 0
		for ioni in range(self.N): 
			for ionj in range(self.N): 
				if ioni != ionj: # skip in case it's the same atom in original unit cell
					dist = np.linalg.norm(self.pos[ioni,] - self.pos[ionj,])
					esum += ( self.get_charges_mult(ioni,ionj) * self.made_const / dist )
		esum *= 14.399645351950543 / 2 # electrostatic constant
		return esum


class Buckingham(Potential):
	'''
	 Calculations for the Buckingham energy contribution.
	'''
	def set_parameters(self, filename, real_cut_off=2):
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

	def get_cutoff(self, hi):
		'''
		 Find how many cells away to check
		 using the minimum cell vector value
		'''
		cutoff = math.ceil(hi/min(self.vects[self.vects!=0]))
		return cutoff

	def calc(self, esum=0):
		'''
		 Interatomic potential
		'''
		chemical_symbols = self.atoms.get_chemical_symbols()
		for ioni in range(self.N):
			for ionj in range(self.N):
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
					if (dist < self.buck[pair]['hi']) & (ioni != ionj):
						esum +=  A*math.exp(-1.0*dist/rho) - C/dist**6
					
					# Check interactions with neighbouring cells
					cutoff = self.get_cutoff(self.buck[pair]['hi'])
					shifts = self.get_shifts( cutoff,self.vects )
					for shift in shifts:
						dist = np.linalg.norm(self.pos[ioni] + shift - self.pos[ionj])
						# Check if distance of ions allows interaction 					
						if (dist < self.buck[pair]['hi']):
							esum += A*math.exp(-1.0*dist/rho) - C/dist**6
		return esum/2