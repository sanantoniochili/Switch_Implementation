import numpy as np
cimport numpy as cnp

cdef class Potential:
	cdef double[:,:] get_shifts(self, int cut_off, double[:,:] vects)
	cdef double[:,:] get_reciprocal_vects(self, double[:,:] vects, double volume)
	cdef int[:,:] map_displacement(self, double[:,:] vects, double[:,:] pos, double[:,:] dx)


cdef class Coulomb(Potential):
	"""Calculations for the Coulomb energy contribution. It
	corresponds to the electrostatic forces exercised among entities.
	
	"""
	cdef double alpha, made_const
	cdef double eself, ereal, erecip
	cdef double[:,:] grad, hessian
	cdef cnp.ndarray chemical_symbols
	cdef int real_cut_off, recip_cut_off
	cdef int[:] charges
	cdef bint param_flag

	cpdef int[:,:] map_displacement(self, double[:,:] vects, double[:,:] pos, double[:,:] dx)
	cpdef set_parameters(self, double alpha, \
		double real_cut_off, double recip_cut_off, \
		cnp.ndarray chemical_symbols, int N, 
		charge_dict, str filename=*)			
	cdef double calc_self(self, int N)
	cdef double calc_real(self, double[:,:] pos, double[:,:] vects, int N) except? -1
	cdef double calc_recip(self, double[:,:] pos, double[:,:] vects, int N) except? -1
	cpdef calc_madelung(self, double[:,:] pos, int N)
	cpdef calc(self, atoms=*, double[:,:] pos_array=*, double[:,:] vects_array=*, int N_=*)
	cdef double[:,:] calc_real_drv(self, double[:,:] pos, double[:,:] vects, int N)
	
	cpdef double[:,:] calc_real_drv2(self, double[:,:] pos, double[:,:] vects, int N)
	
	cdef double[:,:] calc_recip_drv(self, double[:,:] pos, double[:,:] vects, int N)
	cpdef double[:,:] calc_drv(self, atoms=*, double[:,:] pos_array=*, double[:,:] vects_array=*, int N_=*)


cdef class Buckingham(Potential):
	"""Calculations for the Buckingham energy contribution. It
	corresponds to the interatomic forces exercised among entities.
	
	"""
	cdef cnp.ndarray chemical_symbols
	cdef bint param_flag
	cdef double e
	cdef double[:,:] grad
	cdef double[:] radii

	cpdef void set_parameters(self, str filename, 
								cnp.ndarray chemical_symbols, radius_lib=*)
	cpdef int[:] catastrophe_check(self, double[:,:] pos, double fraction, int print_flag=*)
	cdef int get_cutoff(self, double[:,:] vects, float hi)
	cpdef calc(self, atoms=*, double[:,:] pos_array=*, double[:,:] vects_array=*, int N_=*)
	cdef double calc_real(self, double[:,:] pos, double[:,:] vects, int N) except? -1
	cdef double[:,:] calc_drv_(self, double[:,:] pos, double[:,:] vects, int N)
	cpdef double[:,:] calc_drv(self, atoms=*, double[:,:] pos_array=*, double[:,:] vects_array=*, int N_=*)


cdef class Lagrangian(Potential):
	"""Calculations for the Coulomb energy contribution. It
	corresponds to the electrostatic forces exercised among entities.
	
	"""
	cdef double[:,:] llambda, slack
	cdef double econstrain
	cdef double[:,:] grad
	cdef double[:] radii
	cdef bint param_flag

	cpdef set_parameters(self, double[:,:] llambda, str radius_lib, cnp.ndarray chemical_symbols)			
	cpdef double calc_constrain(self, double[:,:] pos, int N) except? -1
	cdef double[:,:] calc_constrain_drv_(self, double[:,:] pos, int N)
	cpdef double[:,:] calc_constrain_drv(self, atoms=*, double[:,:] pos_array=*, int N_=*)