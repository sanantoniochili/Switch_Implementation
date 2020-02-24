# import torch
import sys
import argparse
import fileinput
import numpy as np
from scipy.special import erfc

from cmath import pi
from cmath import exp
import cmath

from ase import *
from ase.io import read as aread
from ase.geometry import Cell

charges = {
	'O': -2.,
	'Sr': 2.,
	'Ti': 4.}

class Coulomb:
	def __init__(self, alpha, real_cut_off=4, recip_cut_off=4):
		self.real_cut_off  = real_cut_off
		self.recip_cut_off = recip_cut_off
		self.alpha         = alpha

	def set_structure(self, atoms, volume, vects, charge_dict):
		'''
		 atoms:       ASE object
		 vects:       Supercell vectors
		 pos:         ASE positions of atoms
		 charge_dict: Dictionary of chemical symbols with charge value
		'''
		self.atoms   = atoms
		self.vects   = vects
		self.pos     = atoms.get_scaled_positions()@vects
		self.volume  = volume
		self.charges = [library[x] for x in atoms.get_chemical_symbols()]

	def get_reciprocal_vects(self):
		'''
		 Calculate reciprocal vectors
		'''
		rvects = np.zeros((3,3))
		for i in np.arange(3):
			rvects[i,] = 2*math.pi*np.cross( self.vects[(1+i)%3,], \
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

		for i in range(0,(2*cut_off+1)**3):
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

	# def calc_self(self, esum):
		'''
		 Calculate self interaction term
		'''
	#     for i in np.arange(N):
	#         esum[i, i] = esum[i, i]-alpha/math.sqrt(math.pi)

	def calc_real(self, vects, esum):
		'''
		 Calculate short range
		'''
		shifts = get_shifts( self.real_cut_off, \
											self.vects )

		for atomi in range(0,N):
			for atomj in range(atomi, N):

				if atomi != atomj: # skip in case it's the same atom
					rij = np.linalg.norm(self.pos[atomi,] - self.pos[atomj,])
					esum[atomi, atomj] += ( get_charges_mult(atomi,atomj) ) * \
												math.erfc( self.alpha*rij )/(2*rij)

				# take care of the rest lattice (+ Ln)
				for shift in shifts:
					rij = np.linalg.norm(self.pos[atomi,] + shift - self.pos[atomj,])
					esum[atomi, atomj] += ( get_charges_mult(atomi,atomj) ) * \
												math.erfc( self.alpha*rij )/(2*rij)

	def calc_recip(self, recip_vects, esum):
		'''
		 Calculate long range
		'''
		shifts = get_shifts( self.recip_cut_off, \
									self.get_reciprocal_vects )
		for atomi in range(0,N):
			for atomj in range(atomi,N):

				dist = self.pos[atomi,] - self.pos[atomj,]

				for s in range(0,shifts):
					k = shifts[s,]

					# -k^2/4a^2
					po = -np.dot(k,k)/(4*alpha**2)
					# 4pi^2 * e^(power) * cos(k(ri-rj))
					numerator = (4*math.pi**2) * math.exp(po) * math.cos(np.dot(k, dist))
					# k^2 * 2piV
					denominator = np.dot(k,k) * (2*math.pi*cell_volume)

					esum[atomi, atomj] += (( get_charges_mult(atomi,atomj) ) * \
														(numerator/denominator))
	
	def calc_complete(self, esum):
		'''
		 Complete calculations
		'''

		# Unit conversion
		esum = esum*14.399645351950543

		# symmetry completion
		for atomi in range(0,N):
			for atomj in range(0,i):
				esum[atomi, atomj] = esum[atomj, atomi]




Atoms = aread("/users/phd/tonyts/Desktop/Data/RandomStart_Sr3Ti3O9/1.cif")


volume = abs(np.linalg.det(vects))
pot = Coulomb(2/(volume**(1.0/3)))

# charges = [library[x] for x in Atoms.get_chemical_symbols()]
# Atoms.set_initial_charges(charges)

# ''' Buckingham-Coulomb '''

# rs = Atoms.get_positions()
# rCell = 2*pi*Atoms.get_reciprocal_cell()  # 2pi not included in ase
# vol = Atoms.get_volume()
# dists = Atoms.get_all_distances()

# '''
#  First version as in 
#  "Ewald Summation for Coulomb Interactions in a Periodic Supercell" p.7

#  - phi = potential field 
#  - e0  = vacuum permitivity
# '''
# e0 = 8.854 * 10**(-12)
# sigma = 0.5
# sqr2 = 2**(1/2)
# sqrp = pi**(1/2)

# '''
#  Short distance energy (Real)
# '''
# Es = 0

# for atomi in range(0, len(dists)):     # distances of dists[atomi] from
#                                        # all other atoms
#     for di in range(0, len(dists[atomi])):
#         if atomi == di:                # exclude distance ri-ri=0
#             continue
#         x = dists[atomi][di] / (sqr2*sigma)
#         Es += 1/(8*pi*e0) * \
#         		((charges[atomi]*charges[di])/dists[atomi][di]) * erfc(x)


# '''
#  Long distance energy (iFF)
# '''
# El = 0                                # long distance energy
#                                       # for all nonzero reciprocal vectors
# for k in rCell:
# 	if ( k==[0,0,0] ).all():
# 		continue
# 	if not ( len(rs)==len(charges) ): # something wrong with atom number
# 		raise ValueError
# 	Sk = 0
# 	for atomi in range(0,len(rs)):     # structure factor Sk
# 		Sk += charges[atomi]*exp( np.dot(k,rs[atomi])*1j )
# 	kn = np.linalg.norm(k)            # get k's norm
# 	Sk2 = (Sk.real**2) + (Sk.imag**2)     # Sk^2
# 	power = - ((sigma**2)/2) * ((kn**2)/2)
# 	El += (1/(2*vol*e0)) * (exp(power)/(kn**2)) * Sk2
# El = El.real
# if not El.imag == 0:
# 	raise ValueError


# '''
#  Self energy (Real)
# '''
# E_ = 0

# for atomi in range(0,len(charges)):
# 	E_ += 1/(4*pi*e0*sqr2*sqrp*sigma) * (charges[atomi]**2)


# '''
#  Result
# '''
# res = Es + El - E_
# print("Es:"+str(Es)+" El: "+str(El)+" Eself: "+str(E_))
# print("Result: "+str(res))
# print(charges)

# https://github.com/SINGROUP/Pysic/blob/master/fortran/Geometry.f90
# https://github.com/vlgusev/IPCSP/blob/master/tools/matrix_generator.py?
