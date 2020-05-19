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
		self.vects = atoms.get_cell()
		self.pos = atoms.get_positions()
		self.volume = abs(np.linalg.det(self.vects))
		self.charges = [charge_dict[x] for x in atoms.get_chemical_symbols()]
		self.N = len(self.pos)
		self.atoms = atoms

	def set_parameters(self):
		pass

	def get_shifts(self, cut_off, vects):
		''' 
		 Get all possible lattice positions:   
		 (2cut_off+1)^3 - {case of (cut_off,cut_off,cut_off)}
		 combinations in R^3  
		'''
		shifts = np.zeros(((2*cut_off+1)**3 - 1, 3))
		tmp = np.array([cut_off, cut_off, cut_off])

		i = 0
		for shift in np.ndindex(2*cut_off+1, 2*cut_off+1, 2*cut_off+1):
			if shift != (cut_off, cut_off, cut_off):
				shifts[i, ] = shift
				shifts[i, ] = shifts[i, ] - tmp
				i = i+1
		shifts = shifts@vects
		return shifts


class Coulomb(Potential):
	'''
	 Calculations for the Coulomb energy contribution.
	 Ewald summation method used for long range.
	'''

	def set_parameters(self, alpha, real_cut_off, \
			recip_cut_off, filename=None):
		"""Set following parameters:

		real_cut_off 	: Images of ions are included up to the 
		larger integer closest this number x real lattice vectors
		per dimension
		recip_cut_off	: Same as real space but number x reciprocal
		lattice vectors
		alpha			: Constant in erfc that controls balance 
		between reciprocal and real space term contribution
		made_const		: Madelung constant if it is to be used
		
		also:
		filename		: Library file for Madelung way potential

		"""
		self.real_cut_off = math.ceil(real_cut_off)
		self.recip_cut_off = math.ceil(recip_cut_off)
		self.alpha = alpha
		self.made_const = 0

		try:
			with open(filename, "r") as fin:
				for line in fin:
					line = line.split()
					found = True
					for symbol in self.atoms.get_chemical_symbols():
						if symbol not in line:
								# chemical formula does not match
							found = False
							break
					if found == True:  # correct formula
						self.made_const = float(line[-1])
		except IOError:
			print("No Madelung library file found.")

	def get_reciprocal_vects(self):
		'''
		 Calculate reciprocal vectors
		'''
		rvects = np.zeros((3, 3))
		for i in np.arange(3):
			rvects[i, ] = 2*pi*np.cross(self.vects[(1+i) % 3, ],
										self.vects[(2+i) % 3]) / self.volume
		return rvects

	def get_charges_mult(self, index1, index2):
		'''
		 Find respective charges and
		 return their product
		'''
		return (self.charges[index1]*self.charges[index2])

	def calc_self(self, esum=[]):
		'''
		 Calculate self interaction term
		'''
		if esum == []:
			esum = np.zeros((self.N, self.N))
		for i in range(0, self.N):
			esum[i, i] -= (self.get_charges_mult(i, i) *
						   (self.alpha / math.sqrt(pi)))
		return esum


	def calc_real_move(self, displace_dict, esum=[]):
		"""Calculate short range when moving single ion

		dx		  : Array that holds the value of displacement
		coord 	  : The coordinate of the ion to be changed
		img_t 	  : Selected image of which an ion will be moved (tuple)
		image 	  : Selected image of which an ion will be moved
		multiplied by the lattice vectors
		ion 	  : Selected ion to be moved
		move_flag : Denotes if an ion was moved

		"""
		dx = np.zeros(3)
		coord = displace_dict['coordinate'] # define affected dimension
		dx[coord]  = displace_dict['value'] # define displacement value
		image = displace_dict['image']
		img_t = tuple(image)
		image = image@self.vects
		ion   = displace_dict['ion']
		move_flag = False

		# energy calculation
		if esum == []:
			esum = np.zeros((self.N, self.N))
		shifts = self.get_shifts(self.real_cut_off, self.vects)
		for ioni in range(0, self.N):
			move_flag = False 
			for ionj in range(ioni, self.N):
				if ioni != ionj:  # skip in case it's the same ion in original unit cell
					# check if it is the correct image
					if (image==[0,0,0]).all():
						# check if it is the correct ion
						if (ioni==ion): # add displacement
							rij = (self.pos[ioni, ]+dx) - self.pos[ionj, ]
							move_flag = True
						elif (ionj==ion):
							rij = self.pos[ioni, ] - (self.pos[ionj, ]+dx)
					else:
						rij = self.pos[ioni, ] - self.pos[ionj, ]
					dist = np.linalg.norm(rij)
					esum[ioni, ionj] += (self.get_charges_mult(ioni, ionj) *
										 math.erfc(self.alpha*dist)/(2*dist))
				# take care of the rest lattice (+ Ln)
				for shift in shifts:
					# check if it is the correct image
					if (ioni != ionj) and (shift==image).all(): # skip same ion
						# check if it is the correct ion
						if (ioni==ion): # add displacement
							rij = (self.pos[ioni, ]+dx) + shift - self.pos[ionj, ]
							move_flag = True
						elif (ionj==ion):
							rij = self.pos[ioni, ] + shift - (self.pos[ionj, ]+dx)
					else:
						rij = self.pos[ioni, ] + shift - self.pos[ionj, ]
					dist = np.linalg.norm(rij)
					esum[ioni, ionj] += (self.get_charges_mult(ioni, ionj) *
										 math.erfc(self.alpha*dist)/(2*dist))

			if move_flag:
				print("Ion {} of image {} moved {}".format(ioni,img_t,dx))

		return esum


	def calc_recip_move(self, recip_vects, displace_dict, esum=[]):
		'''
		 Calculate long range when moving single ion
		'''
		dx = np.zeros(3)
		coord = displace_dict['coordinate'] # define affected dimension
		dx[coord]  = displace_dict['value'] # define displacement value
		image = displace_dict['image']
		img_t = tuple(image)
		image = image@recip_vects
		ion   = displace_dict['ion']

		if esum == []:
			esum = np.zeros((self.N, self.N))
		shifts = self.get_shifts(self.recip_cut_off, recip_vects)
		for ioni in range(0, self.N):
			move_flag = False 
			for ionj in range(ioni, self.N):

				for k in shifts:
					# check if it is the correct image
					if (ioni != ionj) and(k==image).all(): # add displacement
						# check if it is the correct ion
						if (ioni==ion): # add displacement
							rij = (self.pos[ioni, ]+dx) - self.pos[ionj, ]
							move_flag = True
						elif (ionj==ion):
							rij = self.pos[ioni, ] - (self.pos[ionj, ]+dx)
					else:
						rij = self.pos[ioni, ] - self.pos[ionj, ]
					po = -np.dot(k, k)/(4*self.alpha**2)
					numerator = 4 * (pi**2) * (math.exp(po)) * \
						math.cos(np.dot(k, rij))
					denominator = np.dot(k, k) * 2 * pi * self.volume
					esum[ioni, ionj] += ((self.get_charges_mult(ioni, ionj)) *
										 (numerator/denominator))

			if move_flag:
				print("Ion {} of image {} moved {}".format(ioni,img_t,dx))

		return esum

	def calc_complete(self, esum):
		'''
		 Complete lower triangular matrix of monopole to monopole
		'''
		esum *= 14.399645351950543  # electrostatic constant
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
				if ioni != ionj:  # skip in case it's the same atom
												  # in original unit cell
					dist = np.linalg.norm(self.pos[ioni, ] - self.pos[ionj, ])
					esum += (self.get_charges_mult(ioni, ionj)
							 * self.made_const / dist)
		esum *= 14.399645351950543 / 2  # electrostatic constant
		return esum

	def calc_move(self, rvects, displace_dict):
		energies = {}
		Er_array = self.calc_real_move(displace_dict=displace_dict)
		Es_array = self.calc_self()
		Erc_array = self.calc_recip_move(rvects, displace_dict=displace_dict)
		Eupper = Er_array + Es_array + Erc_array

		energies['Real'] = sum(sum(self.calc_complete(Er_array)))
		energies['Self'] = sum(sum(self.calc_complete(Es_array)))
		energies['Reciprocal'] = sum(sum(self.calc_complete(Erc_array)))
		energies['Electrostatic'] = sum(sum(self.calc_complete(Eupper)))
		return energies


class Buckingham(Potential):
	'''
	 Calculations for the Buckingham energy contribution.
	'''

	def set_parameters(self, filename):
		'''
		 Set atom_i-atom_j parameters as read from library file:
		 - par: [A(eV), rho(Angstrom), C(eVAngstrom^6)]
		 - lo : min radius (Angstrom)
		 - hi : max radius (Angstrom)
		'''
		self.buck = {}

		try:
			with open(filename, "r") as fin:
				for line in fin:
					line = line.split()
					if (len(line) < 4):
						continue
					pair = (min(line[0], line[2]), max(line[0], line[2]))
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
		cutoff = math.ceil(hi/min(self.vects[self.vects != 0]))
		return cutoff


	def calc_move(self, displace_dict, esum=0):
		'''
		 Interatomic potential with single ion displacement
		'''
		dx = np.zeros(3)
		coord = displace_dict['coordinate'] # define affected dimension
		dx[coord]  = displace_dict['value'] # define displacement value
		image = displace_dict['image']
		img_t = tuple(image)
		image = image@self.vects
		ion   = displace_dict['ion']


		chemical_symbols = self.atoms.get_chemical_symbols()
		for ioni in range(self.N):
			move_flag = False
			for ionj in range(self.N):
				# Find the pair we are examining
				pair = (min(chemical_symbols[ioni], chemical_symbols[ionj]),
						max(chemical_symbols[ioni], chemical_symbols[ionj]))
				if (pair in self.buck):
					# Pair of ions is listed in parameters file
					A = self.buck[pair]['par'][0]
					rho = self.buck[pair]['par'][1]
					C = self.buck[pair]['par'][2]

					if (ioni != ionj): # skip same ion
						# check if it is the correct image
						if ([0,0,0]==image).all(): # add displacement
							# check if it is the correct ion
							if (ioni==ion): # add displacement
								rij = (self.pos[ioni, ]+dx) - self.pos[ionj, ]
								move_flag = True
							elif (ionj==ion):
								rij = self.pos[ioni, ] - (self.pos[ionj, ]+dx)
						else:
							rij = self.pos[ioni, ] - self.pos[ionj, ]
						dist = np.linalg.norm(rij)
						# Check if distance of ions allows interaction
						if (dist < self.buck[pair]['hi']):
							esum += A*math.exp(-1.0*dist/rho) - C/dist**6

					# Check interactions with neighbouring cells
					cutoff = self.get_cutoff(self.buck[pair]['hi'])
					shifts = self.get_shifts(cutoff, self.vects)
					for shift in shifts:
						# check if it is the correct image
						if (ioni != ionj) and (shift==image).all(): # add displacement
							# check if it is the correct ion
							if (ioni==ion): # add displacement
								rij = (self.pos[ioni, ]+dx) + shift - self.pos[ionj, ]
								move_flag = True
							elif (ionj==ion):
								rij = self.pos[ioni, ] + shift - (self.pos[ionj, ]+dx)
						else:
							rij = self.pos[ioni, ] + shift - self.pos[ionj]
						dist = np.linalg.norm(rij)
						# Check if distance of ions allows interaction
						if (dist < self.buck[pair]['hi']):
							esum += A*math.exp(-1.0*dist/rho) - C/dist**6

			if move_flag:
				print("Ion {} of image {} moved {}".format(ioni,img_t,dx))

		return esum/2
