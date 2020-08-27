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
    def calc_real(self, pos, vects, N):
        """Calculate short range forces
        
        """
        alpha = self.potential.alpha
        a2pi = 2*alpha/pi**(1/2)  # 2a/sqrt(pi)
        grad = np.zeros((N, 3))
        shifts = self.potential.get_shifts(
            self.potential.real_cut_off, vects)

        for ioni in range(0, N):

            ################### PRINTS ####################
            # print("----> Ion: {}{}{} gets the real forces:".format(self.potential.charges[ioni],self.potential.chemical_symbols[ioni],ioni))
            ###############################################

            for ionj in range(0, N):
                if ioni != ionj:  # skip in case it's the same atom or it is constant
                    # direction matters
                    rij = pos[ioni, ] - pos[ionj, ]
                    rnorm = np.linalg.norm(rij)
                    csum = a2pi*math.exp(-alpha**2 * rnorm**2) + \
                        (math.erfc(alpha*rnorm) / rnorm)
                    drv = -self.potential.get_charges_mult(ioni, ionj) *\
                        (rij/(rnorm**2)) * csum  # partial derivative for ion i
                    
                    grad[ioni, ] += -drv
                    grad[ionj, ] -= -drv


                    ################### PRINTS ####################
                    # print("-- from ion {}{}{} of {} eV".format(self.potential.charges[ionj],self.potential.chemical_symbols[ionj],ionj,-self.potential.get_charges_mult(ioni, ionj) *\
                    #             (rij/(rnorm**2)) * csum))
                    ###############################################

                    # take care of the rest lattice (+ Ln)
                    for shift in shifts:
                        rij = pos[ioni, ] + shift - pos[ionj, ]
                        rnorm = np.linalg.norm(rij)
                        csum = a2pi*math.exp(-alpha**2 * rnorm**2) + \
                            (math.erfc(alpha*rnorm) / rnorm)
                        drv = -self.potential.get_charges_mult(ioni, ionj) *\
                            (rij/(rnorm**2)) * csum  # partial derivative for ion i
                        
                        grad[ioni, ] += drv
                        grad[ionj, ] -= drv

                        ################### PRINTS ####################
                        # print("-- in IMAGE {} of {} eV".format(shift,-self.potential.get_charges_mult(ioni, ionj) *\
                        # (rij/(rnorm**2)) * csum))
        # print("----> Total real:")
        # print(-forces/2)
        ###############################################################
        grad = grad/2 * 14.399645351950543  # Coulomb constant
        return grad

    def calc_recip(self, pos, vects, N):
        """Calculate long range forces
        
        """
        alpha = self.potential.alpha
        volume = abs(np.linalg.det(vects))
        grad = np.zeros((N, 3))
        recip_vects = self.potential.get_reciprocal_vects(volume, vects)
        shifts = self.potential.get_shifts(
            self.potential.recip_cut_off, recip_vects)

        for ioni in range(0, N):
            for ionj in range(0, N):
                rij = pos[ioni, ] - pos[ionj, ]
                for k in shifts:
                    po = -np.dot(k, k)/(4*alpha**2)
                    numerator = 4 * (pi**2) * (math.exp(po)) * \
                        k * math.sin(np.dot(k, rij))
                    denominator = np.dot(k, k) * pi * volume
                    drv = -((self.potential.get_charges_mult(ioni, ionj)) *
                                       (numerator/denominator))
                    
                    grad[ioni, ] += drv
                    grad[ionj, ] -= drv


        # print("----> Total recip:")
        # print(-forces/2)
        ###############################################################
        grad = grad/2 * 14.399645351950543  # Coulomb constant
        return grad


class DBuckingham(Forces):
    def calc(self, pos, vects, N):
        """Interatomic forces
        
        """
        chemical_symbols = self.potential.chemical_symbols

        grad = np.zeros((N, 3))
        for ioni in range(N):

            ################### PRINTS ####################
            # print("----> Ion: {}{} gets the interatomic forces:".format(self.potential.chemical_symbols[ioni],ioni))
            ###############################################


            for ionj in range(0, N):
                # Find the pair we are examining
                pair = (min(chemical_symbols[ioni], chemical_symbols[ionj]),
                        max(chemical_symbols[ioni], chemical_symbols[ionj]))
                if (pair in self.potential.buck):
                    # Pair of ions is listed in parameters file
                    A = self.potential.buck[pair]['par'][0]
                    rho = self.potential.buck[pair]['par'][1]
                    C = self.potential.buck[pair]['par'][2]

                    if ioni != ionj:
                        rij = pos[ioni] - pos[ionj]
                        dist = np.linalg.norm(rij)
                        # Check if distance of ions allows interaction
                        if (dist < self.potential.buck[pair]['hi']):
                            csum = -  (A/rho) * \
                                math.exp(-1.0*dist/rho) + 6*C/dist**7
                            drv = (rij/dist) * csum
                            
                            grad[ioni, ] += drv
                            grad[ionj, ] -= drv


                        ################### PRINTS ####################
                        # print("-- from ion {}{} of {} eV".format(self.potential.chemical_symbols[ionj],ionj,-(rij/dist) * csum))
                        ###############################################

                        # Check interactions with neighbouring cells
                        cutoff = self.potential.get_cutoff(
                            vects, self.potential.buck[pair]['hi'])
                        shifts = self.potential.get_shifts(
                            cutoff, vects)
                        for shift in shifts: # this has to be for different i,j 
                                             # or it's constant
                            rij = pos[ioni] + shift - pos[ionj]
                            dist = np.linalg.norm(rij)
                            # Check if distance of ions allows interaction
                            if (dist < self.potential.buck[pair]['hi']):
                                csum = - (A/rho) * \
                                    math.exp(-1.0*dist/rho) + 6*C/dist**7
                                drv = (rij/dist) * csum
                                
                                grad[ioni, ] += drv
                                grad[ionj, ] -= drv

                                ################### PRINTS ####################
                            # print("-- in image {} of {} eV".format(shift,-(rij/dist) * csum))
        # print("----> Total buck:")
        # print(-forces/2)
        ###############################################################
        grad = grad/2
        return grad


if __name__ == "__main__":
    pass
