import numpy as np
cimport numpy as cnp
import warnings
import heapq
import csv

from cython.view cimport array as cvarray
from libc.stdlib cimport malloc, free
from cython.parallel import prange,parallel

from libc cimport bool
from libc.stdio cimport *                                                                
from libc.math cimport *
from libc.float cimport *
from libc.limits cimport *
import cython, shutil

cdef double det3_3(double[:,:] arr):
	"""Returns the determinant of an 3x3 matrix.

	"""
	cdef double det
	det = arr[0,0]*(arr[1,1]*arr[2,2]-arr[1,2]*arr[2,1])- \
			arr[0,1]*(arr[1,0]*arr[2,2]-arr[1,2]*arr[2,0])+ \
			arr[0,2]*(arr[1,0]*arr[2,1]-arr[1,1]*arr[2,0])

	return det

cdef class Potential:

	cdef double[:,:] get_shifts(self, int cut_off, double[:,:] vects):
		"""Returns an array of all possible lattice positions:   
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

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cdef double[:,:] get_reciprocal_vects(self, double[:,:] vects, double volume):
		"""Calculate reciprocal vectors. 
		
		Returns the array of 3D vectors.
		
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
	cdef int[:,:] map_displacement(self, 
		double[:,:] vects, double[:,:] pos, double[:,:] dx):
		"""Checks if given displacements of atoms are illegal with 
		respect to unit cell bounds."""

		if pos==None or (vects.shape[1]!=pos.shape[1]):
			raise ValueError("Positions array is empty.")

		cdef int i, ioni, N = pos.shape[0]
		cdef int dims = pos[0].shape[0]
		cdef double[:,:] rvects
		cdef double* projs
		cdef double pproj, volume = abs(det3_3(vects))
		cdef int[:,:] mask 

		mask = cvarray(shape=(N,3), \
						itemsize=sizeof(int), format="i")

		rvects = self.get_reciprocal_vects(vects, 2*pi)
		with nogil, parallel():
			for ioni in range(N):
				for i in range(3):
					mask[ioni][i] = 1				
					pproj = rvects[i][0]*pos[ioni][0] + \
							rvects[i][1]*pos[ioni][1] + \
							rvects[i][2]*pos[ioni][2]
					if (pproj<0) or (pproj>volume):
						mask[ioni][i] = 0
						printf("Found ion exceeding unit cell bounds\n")
		return mask

cdef class Coulomb(Potential):
	"""Calculations for the Coulomb electrostatic energy contribution.
	The Ewald summation method used for long range. Each sum is calculated
	on a NxN matrix, where N the number of atoms in the unit cell. First 
	the upper triangular matrix is evaluated and the rest is merely repeated,
	thanks to the symmetry of the interactions' effect. Class members:

		real_cut_off 		: Images of ions are included up to the 
		larger integer closest this number x real lattice vectors
		per dimension
		recip_cut_off		: Same as real space but number x reciprocal
		lattice vectors
		alpha				: Constant in erfc that controls balance 
		between reciprocal and real space term contribution
		made_const			: Madelung constant if it is to be used
		charges 			: List of ions' charges in respective positions
		chemical_symbols 	: Ions' chemical symbols in resp. positions 

	"""
	def __init__(self):
		self.alpha = 0
		self.real_cut_off = 0
		self.recip_cut_off = 0
		self.made_const = 0

		self.chemical_symbols = None
		self.charges = None
		self.grad = None
		self.hessian = None

		self.param_flag = False

	cpdef int[:,:] map_displacement(self, 
		double[:,:] vects, double[:,:] pos, double[:,:] dx):
		return Potential.map_displacement(self, vects, pos, dx)

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
		self.grad = cvarray(shape=(N,3), \
								itemsize=sizeof(double), format="d")
		self.hessian = cvarray(shape=(3*N,3*N), \
								itemsize=sizeof(double), format="d")

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
		self.param_flag = True
				
	@cython.boundscheck(False)
	@cython.wraparound(False)
	cdef double calc_self(self, int N): 							
		"""Calculate self interaction term.
		Returns the calculated energy as a float number.

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
		"""Calculate short range energy.
		
		Returns the calculated energy as a float number.
		
		"""
		if pos.shape[1]!=3 or vects.shape[1]!=3:
			raise IndexError("Points are not 3-dimensional.")

		cdef double dist, ereal = 0
		cdef double** esum
		cdef int ioni, ionj, shift, set_excepts
		cdef double[:,:] shifts = self.get_shifts(self.real_cut_off, vects)
		cdef int no_shifts = shifts.shape[0] # number of unit cell images-1
		cdef int value_flag = 0 	

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
					# Avoid division by zero
					if dist == 0:
						dist = DBL_EPSILON
					esum[ioni][ionj] += (self.charges[ioni]*self.charges[ionj] *
											erfc(self.alpha*dist)/(2*dist))
				# take care of the rest lattice (+ Ln)
				for shift in range(no_shifts):
					dist = (pos[ioni, 0]+shifts[shift, 0]-pos[ionj, 0])*(pos[ioni, 0]+shifts[shift, 0]-pos[ionj, 0])+ \
							(pos[ioni, 1]+shifts[shift, 1]-pos[ionj, 1])*(pos[ioni, 1]+shifts[shift, 1]-pos[ionj, 1])+ \
							(pos[ioni, 2]+shifts[shift, 2]-pos[ionj, 2])*(pos[ioni, 2]+shifts[shift, 2]-pos[ionj, 2])
					dist = sqrt(dist)
					# Avoid division by zero
					if dist == 0:
						dist = DBL_EPSILON
					esum[ioni][ionj] += (self.charges[ioni]*self.charges[ionj] *
												erfc(self.alpha*dist)/(2*dist))
		
		# Fill lower triangular matrix with symmetric values
		for ioni in prange(N, nogil=True, schedule='static'):
			for ionj in range(0, ioni):
				esum[ioni][ionj] = esum[ionj][ioni]
			for ionj in range(N):
				ereal += esum[ioni][ionj]

		# Deallocation
		for ioni in prange(N, nogil=True, schedule='static'):
			free(esum[ioni])
		free(esum)

		ereal = ereal*14.399645351950543  # electrostatic constant
		self.ereal = ereal
		return ereal

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cdef double calc_recip(self, double[:,:] pos, double[:,:] vects, int N) except? -1:
		"""Calculate long range energy in reciprocal space.
		
		Returns the calculated energy as a float number.
		
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

						# # avoid underflow
						# if k_2>-log(DBL_EPSILON):
						# 	# term = DBL_EPSILON
						# 	term = 0
						# else:
						term = exp(-k_2/(4*alpha**2))
						# actual calculation
						frac = 4*(pi**2)*term*cos(krij) / (k_2*2*pi*volume)
						esum[ioni][ionj] += self.charges[ioni]*self.charges[ionj]*frac
			free(rij)

		# Fill lower triangular matrix with symmetric values
		for ioni in prange(N, nogil=True, schedule='static'):
			for ionj in range(0, ioni):
				esum[ioni][ionj] = esum[ionj][ioni]
			for ionj in range(N):
				erecip += esum[ioni][ionj]
		
		# Deallocation
		for ioni in prange(N, nogil=True, schedule='static'):
			free(esum[ioni])
		free(esum)

		erecip = erecip*14.399645351950543  # electrostatic constant
		self.erecip = erecip
		return erecip

	@cython.boundscheck(False)
	cpdef calc_madelung(self, double[:,:] pos, int N):
		"""Calculate electrostatic energy using Madelung constant.
		
		Returns the calculated energy as a float number.

		"""
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

	cpdef calc(self, atoms=None, \
		double[:,:] pos_array=None, double[:,:] vects_array=None, int N_=0):
		"""This function needs either the whole Atoms object or
		named arguments for positions (ion positions), vects (unit cell vectors)
		and N (number of atoms in unit cell).
		
		Returns a dictionary with the 3 electrostatic energy parts and their sum.

		"""
		cdef double[:,:] positions 
		cdef double[:,:] vects
		cdef int N

		if atoms:
			positions = atoms.positions
			vects = np.array(atoms.get_cell())
			N = len(positions)
		else:
			positions = pos_array
			vects = vects_array
			N = N_		

		if not self.param_flag:
			raise ValueError("Coulomb potential parameters are not set.")

		self.calc_real(positions,vects,N)
		self.calc_recip(positions,vects,N)
		self.calc_self(N)

		energies = {}
		energies['Real'] = self.ereal
		energies['Self'] = self.eself
		energies['Reciprocal'] = self.erecip
		energies['Electrostatic'] = self.ereal+self.erecip+self.eself

		return energies

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cdef double[:,:] calc_real_drv(self, double[:,:] pos, double[:,:] vects, int N):
		"""Calculate short range electrostatic forces in form of the
		energy function's gradient (forces = -gradient). Instead of 
		using a double N loop, the derivatives of i and j ions are
		updated concurrently.

		Returns the real energy gradient vector.
		
		"""
		if pos.shape[1]!=3 or vects.shape[1]!=3:
			raise IndexError("Points are not 3-dimensional.")

		cdef double dist, dist_2, a2pi, drv, term, alpha = self.alpha
		cdef double* rij
		cdef int ioni, ionj, dim, shift
		cdef double[:,:] shifts = self.get_shifts(self.real_cut_off, vects)
		cdef int no_shifts = shifts.shape[0] # number of unit cell images-1

		a2pi = 2*alpha/sqrt(pi)		

		with nogil, parallel():
			# allocate memory for distance vector
			rij = <double *> malloc(sizeof(double) * 3)
			for ioni in prange(N, schedule='static'):
				for ionj in range(ioni, N):
					if ioni != ionj:  # skip in case it's the same atom or it is constant
						rij[0] = pos[ioni,0]-pos[ionj,0] # distance vector
						rij[1] = pos[ioni,1]-pos[ionj,1]
						rij[2] = pos[ioni,2]-pos[ionj,2]
						dist_2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]
						dist = sqrt(dist_2)
						# # Avoid underflow
						# if alpha*alpha*dist_2>-log(DBL_EPSILON):
						# 	term = DBL_EPSILON
						# 	# term = 0
						# else:
						term = exp(-alpha*alpha*dist_2)
						drv = -self.charges[ioni]*self.charges[ionj] * \
							(
								a2pi*term / dist_2 + \
								erfc(alpha*dist) / (dist_2*dist)
							) # partial derivative without position vector

						# partial deriv with respect to ioni
						self.grad[ioni][0] += drv*rij[0]*14.399645351950543  # Coulomb constant
						self.grad[ioni][1] += drv*rij[1]*14.399645351950543
						self.grad[ioni][2] += drv*rij[2]*14.399645351950543

						# partial deriv with respect to ionj
						self.grad[ionj][0] -= drv*rij[0]*14.399645351950543
						self.grad[ionj][1] -= drv*rij[1]*14.399645351950543
						self.grad[ionj][2] -= drv*rij[2]*14.399645351950543

						# take care of the rest lattice (+ Ln)
						for shift in range(no_shifts):
							rij[0] = pos[ioni,0]+shifts[shift,0]-pos[ionj,0] # distance vector
							rij[1] = pos[ioni,1]+shifts[shift,1]-pos[ionj,1]
							rij[2] = pos[ioni,2]+shifts[shift,2]-pos[ionj,2]
							dist_2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]
							dist = sqrt(dist_2)
							# # Avoid underflow
							# if alpha*alpha*dist_2>-log(DBL_EPSILON):
							# 	term = DBL_EPSILON
							# 	# term = 0
							# else:
							term = exp(-alpha*alpha*dist_2)
							drv = - self.charges[ioni]*self.charges[ionj] * \
								(
									a2pi*term / dist_2 + \
									erfc(alpha*dist) / (dist_2*dist)
								) # partial derivative

							# partial deriv with respect to ioni
							self.grad[ioni][0] += drv*rij[0]*14.399645351950543  # Coulomb constant
							self.grad[ioni][1] += drv*rij[1]*14.399645351950543
							self.grad[ioni][2] += drv*rij[2]*14.399645351950543

							# partial deriv with respect to ionj
							self.grad[ionj][0] -= drv*rij[0]*14.399645351950543
							self.grad[ionj][1] -= drv*rij[1]*14.399645351950543
							self.grad[ionj][2] -= drv*rij[2]*14.399645351950543
			free(rij)
		return self.grad

	@cython.boundscheck(False)
	@cython.wraparound(False)
	@cython.cdivision(True) # for quick modulo division
	cpdef double[:,:] calc_real_drv2(self, double[:,:] pos, double[:,:] vects, int N):
		"""Calculate short range electrostatic forces in form of the
		energy function's Hessian. 

		Returns the real energy Hessian matrix.
		
		"""
		if pos.shape[1]!=3 or vects.shape[1]!=3:
			raise IndexError("Points are not 3-dimensional.")

		# cdef double dist, dist_2, a2pi, drv, term, alpha = self.alpha
		# cdef double* rij
		cdef int i, j, dim, shift, ionl, ionj
		# cdef double[:,:] shifts = self.get_shifts(self.real_cut_off, vects)
		# cdef int no_shifts = shifts.shape[0] # number of unit cell images-1

		# a2pi = 2*alpha/sqrt(pi)	

		i = 1	

		with nogil, parallel():
			# allocate memory for distance vector
			rij = <double *> malloc(sizeof(double) * 3)
			for i in prange(3*N, schedule='static'):
				for j in range(3*N):
					# Operations on coords of same ion
					if <int>floor(i/3)==<int>floor(j/3):
						ionl = <int>floor(i/3)

						for ionj in range(N):
							rij[0] = pos[ionl,0]-pos[ionj,0] # distance vector
							rij[1] = pos[ionl,1]-pos[ionj,1]
							rij[2] = pos[ionl,2]-pos[ionj,2]
							# dist = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2])

							# self.hessian[i][j] += 
							# self.charges[ionl]*self.charges[ionj]/dist



						# Operations on same coords
						if <int>floor(i%3)==<int>floor(j%3):
							pass
						else:
							pass
					else:
						pass
				printf("\n")
		return self.grad

	@cython.boundscheck(False)
	@cython.wraparound(False)
	cdef double[:,:] calc_recip_drv(self, double[:,:] pos, double[:,:] vects, int N):
		"""Calculate long range electrostatic forces in form of the
		energy function's gradient (forces = -gradient). Instead of using 
		a double N loop, the derivatives of i and j ions are 
		updated concurrently.

		Returns the reciprocal energy gradient vector.
		
		"""
		cdef int ioni, ionj, dim, shift
		cdef double drv, k_2, krij, term, alpha = self.alpha
		cdef double volume = abs(det3_3(vects))
		cdef double* rij
		cdef double[:,:] recip_vects = self.get_reciprocal_vects(vects, volume)
		cdef double[:,:] shifts = self.get_shifts(
			self.recip_cut_off, recip_vects)
		cdef int no_shifts = shifts.shape[0] # number of unit cell images-1

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
						# if k_2>-log(DBL_EPSILON):
						# 	term = DBL_EPSILON
						# 	# term = 0
						# else:
						term = exp(-k_2/(4*alpha**2))
						drv = - self.charges[ioni]*self.charges[ionj] * \
										   4*pi*pi*term*sin(krij)/(k_2*pi*volume)
						
						# partial deriv with respect to ioni
						self.grad[ioni][0] += drv*shifts[shift,0]*14.399645351950543  # Coulomb constant
						self.grad[ioni][1] += drv*shifts[shift,1]*14.399645351950543
						self.grad[ioni][2] += drv*shifts[shift,2]*14.399645351950543

						# partial deriv with respect to ionj
						self.grad[ionj][0] -= drv*shifts[shift,0]*14.399645351950543
						self.grad[ionj][1] -= drv*shifts[shift,1]*14.399645351950543
						self.grad[ionj][2] -= drv*shifts[shift,2]*14.399645351950543
			free(rij)
		return self.grad

	cpdef double[:,:] calc_drv(self, atoms=None, \
		double[:,:] pos_array=None, double[:,:] vects_array=None, int N_=0):
		"""Wrapper function to initialise gradient vector and
		call the functions that calculate real and reciprocal parts
		of the partial derivatives. The self term is constant w.r.t
		the ion positions and its derivative is zero.

		Returns the electrostatic energy gradient vector.

		"""
		cdef double[:,:] positions
		cdef double[:,:] vects
		cdef int ioni, dim, N

		if atoms:
			positions = atoms.positions
			vects = np.array(atoms.get_cell())
			N = len(positions)
		else:
			positions = pos_array
			vects = vects_array
			N = N_

		if not self.param_flag:
			raise ValueError("Coulomb potential parameters are not set.")
		
		for ioni in range(N):
			for dim in range(3):
				self.grad[ioni][dim] = 0
		
		self.calc_real_drv(positions,vects,N)
		self.calc_recip_drv(positions,vects,N)
		return self.grad


buck = {}
cdef class Buckingham(Potential):
	"""Calculations for the Buckingham energy contribution. It
	corresponds to the interatomic forces exercised among entities.
	
	"""
	def __init__(self):
		self.chemical_symbols = None
		self.param_flag = False
		self.grad = None
		self.radii = None

	cpdef void set_parameters(self, str filename, 
								cnp.ndarray chemical_symbols, radius_lib=None):
		"""Set atom_i-atom_j parameters as read from library file:
		 - par: [A(eV), rho(Angstrom), C(eVAngstrom^6)]
		 - lo : min radius (Angstrom)
		 - hi : max radius (Angstrom)

		"""
		self.chemical_symbols = chemical_symbols.copy()
		self.grad = cvarray(shape=(len(chemical_symbols),3), \
								itemsize=sizeof(double), format="d")
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

		radius_dict = {}

		if radius_lib:
			self.radii = cvarray(shape=(len(chemical_symbols),), \
					itemsize=sizeof(double), format="d")
			try:
				with open(radius_lib, "r") as fin:
					for line in fin:
						line = line.split()
						radius_dict[line[0]] = float(line[1])
			except IOError:
				print("No library file found for radius values.")
			for s in range(len(chemical_symbols)):
				self.radii[s] = radius_dict[chemical_symbols[s]]
				# print("Set {} for {}".format(radius_dict[chemical_symbols[s]],
				# 	chemical_symbols[s]))

		self.param_flag = True

	cpdef int[:] catastrophe_check(self, double[:,:] pos, double fraction, int print_flag=0):
		"""Check if there are ions in the unit cell that are too close 
		causing Buckingham catastrophe.

		Returns array of pairs if Buckingham catastrophe was discovered, 
		None otherwise.

		"""
		if self.radii == None:
			raise ValueError("Library file for radii not set.")

		cdef int ioni, ionj, N, count, count_true
		cdef double dist, min_thres, min_dist = np.inf
		cdef double keep_thres = 0
		cdef int[:] pairs

		heap = []
		count = 0
		count_true = 0
		symbols = [self.chemical_symbols[x]+str(x) for x in range(len(pos))]

		N = len(pos)
		for ioni in range(N):
			for ionj in range(ioni, N):
				if ioni!=ionj:
					dist = sqrt((pos[ioni, 0]-pos[ionj, 0])*(pos[ioni, 0]-pos[ionj, 0])+ \
								(pos[ioni, 1]-pos[ionj, 1])*(pos[ioni, 1]-pos[ionj, 1])+ \
								(pos[ioni, 2]-pos[ionj, 2])*(pos[ioni, 2]-pos[ionj, 2]))
					elemi = symbols[ioni]
					elemj = symbols[ionj]
					min_thres = self.radii[ioni]+self.radii[ionj]
					if dist < (fraction*min_thres):
						heapq.heappush(heap, [dist, count, ioni, ionj, True])
						print("Detected {},{} too close".format(elemi,elemj))
						count_true += 1
					else:
						heapq.heappush(heap, [dist, count, ioni, ionj, False])
				count += 1

		# Return catastrophe pairs
		count = 0
		if count_true == 0:
			pairs = None
		else:
			pairs = cvarray(shape=(2*count_true,),
					itemsize=sizeof(int), format="i")
		if print_flag==1:
			f = open("dists_temp.csv", 'w')
		while heap:
			heap_elem = heapq.heappop(heap)
			if print_flag==1:
				wr = csv.writer(f, quoting=csv.QUOTE_ALL)
				wr.writerow(heap_elem)
			if heap_elem[4]:
				pairs[count] = heap_elem[2]
				pairs[count+1] = heap_elem[3]
				count += 2
		return pairs


	cdef int get_cutoff(self, double[:,:] vects, float hi):
		"""Find how many cells away to check
		 using the minimum cell vector value (CUBOIDS ONLY)

		 Returns the distance described.
		
		"""
		cdef int cutoff
		cdef double min_cell = INT_MAX
		for x in range(3):
			for y in range(3):
				if (vects[x,y]<min_cell) & (vects[x,y]!=0):
					min_cell = vects[x,y]
		cutoff = int(ceil(hi/min_cell))
		return cutoff

	cpdef calc(self, atoms=None, \
		double[:,:] pos_array=None, double[:,:] vects_array=None, int N_=0):
		"""Interatomic energy potential wrapper.

		Returns the interatomic energy as a float number.
		
		"""
		cdef double[:,:] positions 
		cdef double[:,:] vects
		cdef int N

		if atoms:
			positions = atoms.positions
			vects = np.array(atoms.get_cell())
			N = len(positions)
		else:
			positions = pos_array
			vects = vects_array
			N = N_	

		if not self.param_flag:
			raise ValueError("Buckingham potential parameters are not set.")

		return self.calc_real(positions, vects, N)

	cdef double calc_real(self, double[:,:] pos, double[:,:] vects, int N) except? -1:
		"""Fucnction to calculate the Buckingham potential.

		Returns the interatomic energy as a float number.

		"""
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

	cdef double[:,:] calc_drv_(self, double[:,:] pos, double[:,:] vects, int N):
		"""Interatomic forces in form of the
		energy function's gradient (forces = -gradient). 
		Instead of using a double N loop,
		the derivatives of i and j ions are updated concurrently.

		Returns the interatomic energy gradient vector.
		
		"""
		if not self.param_flag:
			raise ValueError("Buckingham potential parameters are not set.")

		cdef double* rij = <double *> malloc(sizeof(double) * 3)     
		for ioni in range(N):
			for ionj in range(ioni, N):
				# Find the pair we are examining
				pair = (min(self.chemical_symbols[ioni], self.chemical_symbols[ionj]),
						max(self.chemical_symbols[ioni], self.chemical_symbols[ionj]))
				if (pair in buck):
					# Pair of ions is listed in parameters file
					A = buck[pair]['par'][0]
					rho = buck[pair]['par'][1]
					C = buck[pair]['par'][2]

					if ioni != ionj:
						rij[0] = pos[ioni,0]-pos[ionj,0] # distance vector
						rij[1] = pos[ioni,1]-pos[ionj,1]
						rij[2] = pos[ioni,2]-pos[ionj,2]
						dist = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2])
						# Check if distance of ions allows interaction
						if (dist < buck[pair]['hi']):
							drv = -  (A/rho) * \
								exp(-1.0*dist/rho) + 6*C/dist**7
							
							# partial deriv with respect to ioni
							self.grad[ioni][0] += drv*rij[0]/dist
							self.grad[ioni][1] += drv*rij[1]/dist
							self.grad[ioni][2] += drv*rij[2]/dist

							# partial deriv with respect to ionj
							self.grad[ionj][0] -= drv*rij[0]/dist
							self.grad[ionj][1] -= drv*rij[1]/dist
							self.grad[ionj][2] -= drv*rij[2]/dist

						# Check interactions with neighbouring cells
						cutoff = self.get_cutoff(vects, buck[pair]['hi'])
						shifts = self.get_shifts(cutoff, vects)
						no_shifts = shifts.shape[0] # number of unit cell images-1
						for shift in range(no_shifts):
							rij[0] = pos[ioni,0]+shifts[shift,0]-pos[ionj,0] # distance vector
							rij[1] = pos[ioni,1]+shifts[shift,1]-pos[ionj,1]
							rij[2] = pos[ioni,2]+shifts[shift,2]-pos[ionj,2]
							dist = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2])
							# Check if distance of ions allows interaction
							if (dist < buck[pair]['hi']):
								drv = - (A/rho) * \
									exp(-1.0*dist/rho) + 6*C/dist**7
								
								# partial deriv with respect to ioni
								self.grad[ioni][0] += drv*rij[0]/dist
								self.grad[ioni][1] += drv*rij[1]/dist
								self.grad[ioni][2] += drv*rij[2]/dist

								# partial deriv with respect to ionj
								self.grad[ionj][0] -= drv*rij[0]/dist
								self.grad[ionj][1] -= drv*rij[1]/dist
								self.grad[ionj][2] -= drv*rij[2]/dist
		free(rij)
		return self.grad

	cpdef double[:,:] calc_drv(self, atoms=None, \
		double[:,:] pos_array=None, double[:,:] vects_array=None, int N_=0):
		"""Wrapper function to initialise gradient vector and
		call the function that calculates the gradient.

		Returns the interatomic energy gradient vector.

		"""
		cdef double[:,:] positions
		cdef double[:,:] vects
		cdef int ioni, dim, N

		if atoms:
			positions = atoms.positions
			vects = np.array(atoms.get_cell())
			N = len(positions)
		else:
			positions = pos_array
			vects = vects_array
			N = N_	

		if not self.param_flag:
			raise ValueError("Coulomb potential parameters are not set.")
		
		for ioni in range(N):
			for dim in range(3):
				self.grad[ioni][dim] = 0
		
		self.calc_drv_(positions,vects,N)
		return self.grad

cdef class Lagrangian(Potential):
	def __init__(self):
		self.param_flag = False
		self.grad = None

	cpdef set_parameters(self, double[:,:] llambda, str radius_lib, 
		cnp.ndarray chemical_symbols):

		cdef int s, ioni, ionj
		
		# Make sure lambda is a symmetric matrix
		for ioni in range(llambda.shape[0]):
			for ionj in range(ioni, llambda.shape[0]):
				assert(llambda[ioni][ionj]==llambda[ionj][ioni])

		self.llambda = llambda
		self.econstrain = 0
		self.radii = cvarray(shape=(len(chemical_symbols),), \
					itemsize=sizeof(double), format="d")
		radius_dict = {}

		try:
			with open(radius_lib, "r") as fin:
				for line in fin:
					line = line.split()
					radius_dict[line[0]] = float(line[1])
		except IOError:
			print("No library file found for radius values.")

		for s in range(len(chemical_symbols)):
			self.radii[s] = radius_dict[chemical_symbols[s]]

		self.grad = cvarray(shape=(len(chemical_symbols),3), \
								itemsize=sizeof(double), format="d")
		self.param_flag = True
			
	cpdef double calc_constrain(self, double[:,:] pos, int N) except? -1:

		if not self.param_flag:
			raise ValueError("Constraint parameters are not set.")

		self.econstrain = 0

		cdef double dist, ctemp
		cdef int ioni, ionj

		for ioni in range(N):
			for ionj in range(ioni+1, N):
				dist = (pos[ioni, 0]-pos[ionj, 0])*(pos[ioni, 0]-pos[ionj, 0])+ \
						(pos[ioni, 1]-pos[ionj, 1])*(pos[ioni, 1]-pos[ionj, 1])+ \
						(pos[ioni, 2]-pos[ionj, 2])*(pos[ioni, 2]-pos[ionj, 2])
				dist = sqrt(dist)


				# Add calculated constraint to final energy
				# When radii sum is greater it adds cost to minimisation
				ctemp = self.llambda[ioni][ionj]* \
					(self.radii[ioni]+self.radii[ionj]-dist)
				self.econstrain += ctemp

				# # Print deed
				# print("Constraint: \t{} Ions: {},{} \
				# 	\t Distance: {} \t Radii: {},{} ".format(
				# 		ctemp, ioni, ionj, dist, 
				# 		self.radii[iona], self.radii[ionb]))
				# columns = shutil.get_terminal_size().columns
				# for c in range(columns):
				# 	print("-", end="")
				# print()

		return self.econstrain

	cdef double[:,:] calc_constrain_drv_(self, double[:,:] pos, int N):
		if not self.param_flag:
			raise ValueError("Constraint parameters are not set.")

		cdef double dist
		cdef int ioni, ionj
		cdef double* rij = <double *> malloc(sizeof(double) * 3)     
		
		for ioni in range(N):
			for ionj in range(ioni+1, N):

				rij[0] = pos[ioni,0]-pos[ionj,0] # distance vector
				rij[1] = pos[ioni,1]-pos[ionj,1]
				rij[2] = pos[ioni,2]-pos[ionj,2]
				dist = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2])

				# partial deriv with respect to ioni
				self.grad[ioni][0] -= self.llambda[ioni][ionj]*rij[0]/dist
				self.grad[ioni][1] -= self.llambda[ioni][ionj]*rij[1]/dist
				self.grad[ioni][2] -= self.llambda[ioni][ionj]*rij[2]/dist

				# partial deriv with respect to ionj
				self.grad[ionj][0] += self.llambda[ioni][ionj]*rij[0]/dist
				self.grad[ionj][1] += self.llambda[ioni][ionj]*rij[1]/dist
				self.grad[ionj][2] += self.llambda[ioni][ionj]*rij[2]/dist
		free(rij)
		return self.grad

	cpdef double[:,:] calc_constrain_drv(self, atoms=None, \
		double[:,:] pos_array=None, int N_=0):
		"""Wrapper function to initialise gradient 
		vector of the added constraints.

		"""
		cdef double[:,:] positions
		cdef double[:,:] vects
		cdef int ioni, dim, N

		if atoms:
			positions = atoms.positions
			N = len(positions)
		else:
			positions = pos_array
			N = N_

		if not self.param_flag:
			raise ValueError("Lagrange potential parameters are not set.")
		
		for ioni in range(N):
			for dim in range(3):
				self.grad[ioni][dim] = 0
		
		self.calc_constrain_drv_(positions,N)
		return self.grad


