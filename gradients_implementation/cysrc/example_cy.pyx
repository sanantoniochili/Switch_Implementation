import numpy as np
cimport numpy as cnp

from cython.view cimport array as cvarray
from libc.stdlib cimport malloc, free
from cython.parallel import prange,parallel
from libc.math cimport ceil
from libc.stdio cimport *                                                                
from libc.math cimport *
from libc.float cimport *
from libc.limits cimport *
import cython


cdef double det3_3(double[:,:] arr):
	"""Get the determinant of an 3x3 matrix.

	"""
	cdef double det
	det = arr[0,0]*(arr[1,1]*arr[2,2]-arr[1,2]*arr[2,1])- \
			arr[0,1]*(arr[1,0]*arr[2,2]-arr[1,2]*arr[2,0])+ \
			arr[0,2]*(arr[1,0]*arr[2,1]-arr[1,1]*arr[2,0])

	return det

cdef class Potential:

	cdef double[:,:] get_shifts(self, int cut_off, double[:,:] vects):
		"""Get all possible lattice positions:   
		 (2cut_off+1)^3 - {case of (cut_off,cut_off,cut_off)}
		 combinations in R^3  
		
		"""
		cdef int dim = (2*cut_off+1)**3 - 1
		shifts = cvarray(shape=(dim, 3), itemsize=sizeof(double), format="d")
		cdef double [:,:] shifts_view = shifts
		cdef int i = 0


		for no1 in range(2*cut_off+1):
			for no2 in range(2*cut_off+1):
				for no3 in range(2*cut_off+1):
					if (no1!=cut_off) or (no2!=cut_off) or (no3!=cut_off):
						shifts_view[i, 0] = no1
						shifts_view[i, 1] = no2
						shifts_view[i, 2] = no3
						for item in range(3):
							shifts_view[i, item] = shifts_view[i, item] - cut_off
						i = i+1
		shifts_view = np.dot(shifts_view,vects)
		return shifts_view

cdef class Coulomb(Potential):
	"""Calculations for the Coulomb electrostatic energy contribution.
	The Ewald summation method used for long range. Each sum is calculated
	on a NxN matrix, where N the number of atoms in the unit cell. First 
	the upper triangular matrix is evaluated and the rest is merely repeated,
	thanks to the symmetry of the interactions' effect. Class members:

		real_cut_off 	: Images of ions are included up to the 
		larger integer closest this number x real lattice vectors
 		per dimension
 		recip_cut_off	: Same as real space but number x reciprocal
 		lattice vectors
 		alpha			: Constant in erfc that controls balance 
 		between reciprocal and real space term contribution
 		made_const		: Madelung constant if it is to be used
 		charges 		: List of ions' charges in respective positions
 		chemical_symbols: Ions' chemical symbols in resp. positions 

	"""
	cdef double alpha, made_const
	cdef double eself, ereal, erecip
	cdef cnp.ndarray chemical_symbols
	cdef int real_cut_off, recip_cut_off
	cdef int[:] charges

	def __init__(self):
		self.alpha = 0
		self.real_cut_off = 0
		self.recip_cut_off = 0
		self.made_const = 0

		self.chemical_symbols = None
		self.charges = None


	@cython.boundscheck(False)
	@cython.wraparound(False)
	cpdef set_parameters(self, double alpha, \
		double real_cut_off, double recip_cut_off, \
		cnp.ndarray chemical_symbols, int N, 
		charge_dict, str filename=None):
		"""Set the class's parameters and specify the following:
 		N 				: Number of ions
 		filename		: Library file for Madelung way potential

 		"""
		self.real_cut_off = int(ceil(real_cut_off))
		self.recip_cut_off = int(ceil(recip_cut_off))
		self.alpha = alpha
		self.made_const = 0
		self.chemical_symbols = chemical_symbols.copy()
		self.eself = 0
		self.ereal = 0
		self.erecip = 0

		cdef int count = 0
		self.charges = cvarray(shape=(N,), \
								itemsize=sizeof(int), format="i")
		for elem in chemical_symbols: # for every ion get its ch.symbol "elem"
			self.charges[count] = charge_dict[elem]
			count += 1
        
		if filename:
			try:
				with open(filename,"r") as fin:
					for line in fin:
						line = line.split()
						found = True
						for symbol in chemical_symbols:
							if symbol not in line: # chemical formula does not match
								found = False
								break
							if found == True:  # correct formula
								self.made_const = float(line.shape[0]-2)
			except IOError:
				print("No Madelung library file found.")
                
	@cython.boundscheck(False)
	@cython.wraparound(False)
	cdef double[:,:] get_reciprocal_vects(self, double[:,:] vects, double volume):
		"""Calculate reciprocal vectors.
		
		"""
		cdef int i,a,b
		rvects = cvarray(shape=(3,3), \
						itemsize=sizeof(double), format="d")
		for i in range(3):
			a = (1+i) % 3
			b = (2+i) % 3
			rvects[i, ] = 2*pi*np.cross(vects[a, ],
										vects[b, ]) / volume
		return rvects

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cdef double calc_self(self, int N): 							
		"""Calculate self interaction term
		
		"""
		cdef int i
		cdef double eself = 0

		self.eself= 0
		for i in range(N):
			eself -= (self.charges[i]*self.charges[i] *
						   (self.alpha / sqrt(pi)))
		self.eself = eself*14.399645351950543  # electrostatic constant
		return eself

	@cython.boundscheck(False)
	@cython.wraparound(False)	
	cdef double calc_real(self, double[:,:] pos, double[:,:] vects, int N) except? -1:
		"""Calculate short range
		
		"""
		if pos.shape[1]!=3 or vects.shape[1]!=3:
			raise IndexError("Points are not 3-dimensional.")

		cdef double dist, ereal = 0
		cdef double** esum
		cdef int ioni, ionj, shift
		cdef double[:,:] shifts = self.get_shifts(self.real_cut_off, vects)
		cdef int no_shifts = shifts.shape[0] # number of unit cell images-1		

		self.ereal = 0

		# create array with sums for each N*N position
		esum = <double **> malloc(sizeof(double *) * N)
		for ioni in range(N):
			esum[ioni] = <double *> malloc(sizeof(double) * N)
			for ionj in range(N):
				esum[ioni][ionj] = 0

		for ioni in prange(N, nogil=True, schedule='static'):
			for ionj in range(ioni, N):
				if ioni != ionj:  # skip in case it's the same ion in original unit cell
					dist = (pos[ioni, 0]-pos[ionj, 0])*(pos[ioni, 0]-pos[ionj, 0])+ \
							(pos[ioni, 1]-pos[ionj, 1])*(pos[ioni, 1]-pos[ionj, 1])+ \
							(pos[ioni, 2]-pos[ionj, 2])*(pos[ioni, 2]-pos[ionj, 2])
					dist = sqrt(dist)
					esum[ioni][ionj] += (self.charges[ioni]*self.charges[ionj] *
											erfc(self.alpha*dist)/(2*dist))
				# take care of the rest lattice (+ Ln)
				for shift in range(no_shifts):
					dist = (pos[ioni, 0]+shifts[shift, 0]-pos[ionj, 0])*(pos[ioni, 0]+shifts[shift, 0]-pos[ionj, 0])+ \
							(pos[ioni, 1]+shifts[shift, 1]-pos[ionj, 1])*(pos[ioni, 1]+shifts[shift, 1]-pos[ionj, 1])+ \
							(pos[ioni, 2]+shifts[shift, 2]-pos[ionj, 2])*(pos[ioni, 2]+shifts[shift, 2]-pos[ionj, 2])
					dist = sqrt(dist)
					esum[ioni][ionj] += (self.charges[ioni]*self.charges[ionj] *
										 		erfc(self.alpha*dist)/(2*dist))
		# Deallocation
		for ioni in prange(N, nogil=True, schedule='static'):
			for ionj in range(0, ioni):
				esum[ioni][ionj] = esum[ionj][ioni] # complete lower triangular
			for ionj in range(N):
				ereal += esum[ioni][ionj]
			free(esum[ioni])
		free(esum)

		ereal = ereal*14.399645351950543  # electrostatic constant
		self.ereal = ereal
		return ereal

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cdef double calc_recip(self, double[:,:] pos, double[:,:] vects, int N) except? -1:
		"""Calculate long range
		
		"""
		if pos.shape[1]!=3 or vects.shape[1]!=3:
			raise IndexError("Points are not 3-dimensional.")
		
		cdef double* rij 
		cdef double** esum
		cdef int ioni, ionj, shift
		cdef double volume = abs(det3_3(vects))
		cdef double[:,:] rvects = self.get_reciprocal_vects(vects, volume)
		cdef double[:,:] shifts = self.get_shifts(self.recip_cut_off, rvects)
		cdef int no_shifts = shifts.shape[0] # number of unit cell images-1
		cdef double k_2, krij, frac, term, alpha = self.alpha
		cdef double erecip = 0

		self.erecip = 0

		# create array with sums for each N*N position
		esum = <double **> malloc(sizeof(double *) * N)
		for ioni in range(N):
			esum[ioni] = <double *> malloc(sizeof(double) * N)
			for ionj in range(N):
				esum[ioni][ionj] = 0

		with nogil, parallel():
			# allocate memory for distance vector
			rij = <double *> malloc(sizeof(double) * 3)
			for ioni in prange(N, schedule='static'):
				for ionj in range(ioni, N):
					rij[0] = pos[ioni,0]-pos[ionj,0] # distance vector
					rij[1] = pos[ioni,1]-pos[ionj,1]
					rij[2] = pos[ioni,2]-pos[ionj,2]
					for shift in range(no_shifts):
						# shift on 2nd power
						k_2 = shifts[shift, 0]*shifts[shift, 0]+ \
								shifts[shift, 1]*shifts[shift, 1]+ \
								shifts[shift, 2]*shifts[shift, 2]

						# dot product to find image
						krij = shifts[shift, 0]*rij[0]+ \
								shifts[shift, 1]*rij[1]+ \
								shifts[shift, 2]*rij[2]

						# avoid underflow
						if k_2>-log(DBL_EPSILON):
							term = DBL_EPSILON
						else:
							term = exp(-k_2/(4*alpha**2))
						# actual calculation
						frac = 4*(pi**2)*term*cos(krij) / (k_2*2*pi*volume)
						esum[ioni][ionj] += self.charges[ioni]*self.charges[ionj]*frac
			free(rij)

		# Deallocation
		for ioni in prange(N, nogil=True, schedule='static'):
			for ionj in range(0, ioni):
				esum[ioni][ionj] = esum[ionj][ioni] # complete lower triangular
			for ionj in range(N):
				erecip += esum[ioni][ionj]
			free(esum[ioni])
		free(esum)

		erecip = erecip*14.399645351950543  # electrostatic constant
		self.erecip = erecip
		return erecip

	@cython.boundscheck(False)
	cpdef calc_madelung(self, double[:,:] pos, int N):
		if not self.made_const:
			return None
		cdef double dist, esum = 0
		cdef int ioni, ionj

		for ioni in prange(N, nogil=True, schedule='static'):
			for ionj in range(N):
				if ioni != ionj:  # skip in case it's the same atom
												  # in original unit cell
					dist = (pos[ioni, 0]-pos[ionj, 0])*(pos[ioni, 0]-pos[ionj, 0])+ \
							(pos[ioni, 1]-pos[ionj, 1])*(pos[ioni, 1]-pos[ionj, 1])+ \
							(pos[ioni, 2]-pos[ionj, 2])*(pos[ioni, 2]-pos[ionj, 2])
					esum += (self.charges[ioni]*self.charges[ionj]* \
								self.made_const / dist)
		esum *= 14.399645351950543 / 2  # Coulomb constant
		return esum

	cpdef calc(self, atoms):
		"""This function needs either the whole Atoms object or
		named arguments for positions (ion positions), vects (unit cell vectors)
		and N (number of atoms in unit cell)

		"""
		positions = atoms.positions
		vects = np.array(atoms.get_cell())
		N = len(positions)

		self.calc_real(positions,vects,N)
		self.calc_recip(positions,vects,N)
		self.calc_self(N)

		energies = {}
		energies['Real'] = self.ereal
		energies['Self'] = self.eself
		energies['Reciprocal'] = self.erecip
		energies['Electrostatic'] = self.ereal+self.erecip+self.eself
		return energies

buck = {}
cdef class Buckingham(Potential):
	"""Calculations for the Buckingham energy contribution. It
	corresponds to the interatomic forces exercised among entities.
	
	"""
	cdef cnp.ndarray chemical_symbols
	cdef double e

	def __init__(self):
		self.chemical_symbols = None

	cpdef void set_parameters(self, str filename, 
								cnp.ndarray chemical_symbols):
		"""Set atom_i-atom_j parameters as read from library file:
		 - par: [A(eV), rho(Angstrom), C(eVAngstrom^6)]
		 - lo : min radius (Angstrom)
		 - hi : max radius (Angstrom)

		"""
		self.chemical_symbols = chemical_symbols.copy()
		try:
			with open(filename, "r") as fin:
				for line in fin:
					line = line.split()
					if (len(line) < 4):
						continue
					pair = (min(line[0], line[2]), max(line[0], line[2]))
					buck[pair] = {}
					buck[pair]['par'] = list(map(float, line[4:7]))
					buck[pair]['lo'] = float(line[7])
					buck[pair]['hi'] = float(line[-1])
		except IOError:
			print("No library file found for Buckingham constants.")

	cdef int get_cutoff(self, double[:,:] vects, float hi):
		"""Find how many cells away to check
		 using the minimum cell vector value
		
		"""
		cdef int cutoff
		cdef double min_cell = INT_MAX
		for x in range(3):
			for y in range(3):
				if (vects[x,y]<min_cell) & (vects[x,y]!=0):
					min_cell = vects[x,y]
		cutoff = int(ceil(hi/min_cell))
		return cutoff

	cpdef calc(self, atoms):
		"""Interatomic potential
		
		"""
		pos = atoms.get_positions()
		N = len(atoms.get_positions())
		vects = np.array(atoms.get_cell())
		return self.calc_real(pos, vects, N)

	cdef double calc_real(self, double[:,:] pos, double[:,:] vects, int N) except? -1:
		cdef double A,rho, C, dist, esum = 0
		cdef int ioni, ionj, no_shifts

		self.e = 0

		for ioni in range(N):
			for ionj in range(N):
				# Find the pair we are examining
				pair = (min(self.chemical_symbols[ioni], self.chemical_symbols[ionj]),
						max(self.chemical_symbols[ioni], self.chemical_symbols[ionj]))
				if (pair in buck):
					# Pair of ions is listed in parameters file
					A = buck[pair]['par'][0]
					rho = buck[pair]['par'][1]
					C = buck[pair]['par'][2]

					dist = (pos[ioni, 0]-pos[ionj, 0])*(pos[ioni, 0]-pos[ionj, 0])+ \
							(pos[ioni, 1]-pos[ionj, 1])*(pos[ioni, 1]-pos[ionj, 1])+ \
							(pos[ioni, 2]-pos[ionj, 2])*(pos[ioni, 2]-pos[ionj, 2])
					dist = sqrt(dist)
					# Check if distance of ions allows interaction
					if (dist < buck[pair]['hi']) & (ioni != ionj):
						esum += A*exp(-1.0*dist/rho) - C/dist**6

					# Check interactions with neighbouring cells
					cutoff = self.get_cutoff(vects, buck[pair]['hi'])
					shifts = self.get_shifts(cutoff, vects)
					no_shifts = shifts.shape[0] # number of unit cell images-1
					for shift in range(no_shifts):
						dist = (pos[ioni, 0]+shifts[shift, 0]-pos[ionj, 0])*(pos[ioni, 0]+shifts[shift, 0]-pos[ionj, 0])+ \
								(pos[ioni, 1]+shifts[shift, 1]-pos[ionj, 1])*(pos[ioni, 1]+shifts[shift, 1]-pos[ionj, 1])+ \
								(pos[ioni, 2]+shifts[shift, 2]-pos[ionj, 2])*(pos[ioni, 2]+shifts[shift, 2]-pos[ionj, 2])
						dist = sqrt(dist)
						# Check if distance of ions allows interaction
						if (dist < buck[pair]['hi']):
							esum += A*exp(-1.0*dist/rho) - C/dist**6
		self.e = esum/2
		return esum/2