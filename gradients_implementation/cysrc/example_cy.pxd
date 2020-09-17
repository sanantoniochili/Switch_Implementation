# from cmath import pi
from libc.math cimport pi


cdef class Potential:
	cdef double[:,:] get_shifts(self, int cut_off, double[:,:] vects)
cdef class Coulomb:
	cdef double[:,:] get_shifts(self, int cut_off, double[:,:] vects)
	cpdef void set_parameters(self, double alpha, double real_cut_off, double recip_cut_off, \
		numpy.ndarray chemical_symbols, charge_dict, str filename)
	cdef double[:,:] get_reciprocal_vects(self, double[:,:] vects,  double volume)
	cdef void calc_self(self, N)
	cdef int calc_real(self, double[:,:] pos, double[:,:] vects, int N)
	cdef int calc_recip(self, double[:,:] pos, double[:,:] vects, int N)
	cpdef calc_madelung(self, double[:,:] pos, int N)
	def calc(self, atoms)

cdef class Buckingham:
	cpdef void set_parameters(self, str filename, cnp.ndarray chemical_symbols)
	cpdef int get_cutoff(self, double[:,:] vects, float hi)
	cpdef calc(self, atoms)
	cdef double calc_real(self, double[:,:] pos, double[:,:] vects, int N)
