from scipy.special import erfc
import numpy as np

from cmath import pi
from cmath import exp
import cmath
import math

from ase import *
from ase.visualize import view
from ase.geometry import Cell


class Forces:
    def __init__(self, potential):
        self.potential = potential


class DCoulomb(Forces):
    def calc_real(self):
        '''
         Calculate short range forces
        '''
        alpha = self.potential.alpha
        N = self.potential.N
        a2pi = 2*alpha/pi**(1/2)  # 2a/sqrt(pi)

        forces = np.zeros((N, 3))
        shifts = self.potential.get_shifts(
            self.potential.real_cut_off, self.potential.vects)
        for ioni in range(0, N):
            for ionj in range(0, N):

                if ioni != ionj:  # skip in case it's the same atom in original unit cell
                    # direction matters
                    rij = self.potential.pos[ioni, ] - \
                        self.potential.pos[ionj, ]
                    rnorm = np.linalg.norm(rij)
                    csum = a2pi*math.exp(-alpha**2 * rnorm**2) + \
                        (math.erfc(alpha*rnorm) / rnorm)
                    forces[ioni, ] -= (1/2) * self.potential.get_charges_mult(ioni, ionj) *\
                        (rij/(rnorm**2)) * csum  # partial derivative for ion i
                # take care of the rest lattice (+ Ln)
                for shift in shifts:
                    rij = self.potential.pos[ioni, ] + \
                        shift - self.potential.pos[ionj, ]
                    rnorm = np.linalg.norm(rij)
                    csum = a2pi*math.exp(-alpha**2 * rnorm**2) + \
                        (math.erfc(alpha*rnorm) / rnorm)
                    forces[ioni, ] -= (1/2) * self.potential.get_charges_mult(ioni, ionj) *\
                        (rij/(rnorm**2)) * csum  # partial derivative for ion i
        return forces

    def calc_recip(self):
        '''
         Calculate long range forces
        '''
        alpha = self.potential.alpha
        N = self.potential.N

        forces = np.zeros((N, 3))
        recip_vects = self.potential.get_reciprocal_vects()
        shifts = self.potential.get_shifts(
            self.potential.recip_cut_off, recip_vects)
        for ioni in range(0, N):
            for ionj in range(0, N):

                rij = self.potential.pos[ioni, ] - self.potential.pos[ionj, ]
                for k in shifts:
                    po = -np.dot(k, k)/(4*alpha**2)
                    numerator = 4 * (pi**2) * (math.exp(po)) * \
                        k * math.sin(np.dot(k, rij))
                    denominator = np.dot(k, k) * 2 * pi * self.potential.volume
                    forces[ioni, ] -= ((self.potential.get_charges_mult(ioni, ionj)) *
                                       (numerator/denominator))
        return forces


class DBuckingham(Forces):
    def calc(self):
        '''
         Interatomic forces
        '''
        chemical_symbols = self.potential.atoms.get_chemical_symbols()
        N = self.potential.N

        forces = np.zeros((N, 3))
        for ioni in range(N):
            for ionj in range(N):
                # Find the pair we are examining
                pair = (min(chemical_symbols[ioni], chemical_symbols[ionj]),
                        max(chemical_symbols[ioni], chemical_symbols[ionj]))
                if (pair in self.potential.buck):
                    # Pair of ions is listed in parameters file
                    A = self.potential.buck[pair]['par'][0]
                    rho = self.potential.buck[pair]['par'][1]
                    C = self.potential.buck[pair]['par'][2]

                    rij = self.potential.pos[ioni] - self.potential.pos[ionj]
                    dist = np.linalg.norm(rij)
                    # Check if distance of ions allows interaction
                    if (dist < self.potential.buck[pair]['hi']) & (ioni != ionj):
                        csum = -  (A/rho) * \
                            math.exp(-1.0*dist/rho) + 6*C/dist**7
                        forces[ioni] += (rij/dist) * csum
                    # Check interactions with neighbouring cells
                    cutoff = self.potential.get_cutoff(
                        self.potential.buck[pair]['hi'])
                    shifts = self.potential.get_shifts(
                        cutoff, self.potential.vects)
                    for shift in shifts:
                        rij = self.potential.pos[ioni] + \
                            shift - self.potential.pos[ionj]
                        dist = np.linalg.norm(rij)
                        # Check if distance of ions allows interaction
                        if (dist < self.potential.buck[pair]['hi']):
                            csum = - (A/rho) * \
                                math.exp(-1.0*dist/rho) + 6*C/dist**7
                            forces[ioni] += (rij/dist) * csum

        return forces


if __name__ == "__main__":
    pass
