from libc.math cimport pi


cdef class Potential
cdef class Coulomb:
	cpdef void set_parameters(self, double alpha, double real_cut_off, double recip_cut_off, \
		numpy.ndarray chemical_symbols, charge_dict, str filename)
	cpdef calc_madelung(self, double[:,:] pos, int N)
	cpdef calc(self, atoms)
	cpdef double[:,:] calc_real_drv(self, double[:,:] pos, double[:,:] vects, int N)

cdef class Buckingham:
	cpdef void set_parameters(self, str filename, cnp.ndarray chemical_symbols)
	cpdef int get_cutoff(self, double[:,:] vects, float hi)
	cpdef calc(self, atoms)
	cpdef double[:,:] calc_drv(self, atoms)
