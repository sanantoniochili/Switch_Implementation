# cy = timeit.timeit('''example_cy.test(5)''',setup='import example_cy',number=100)
# pi = timeit.timeit('''example_cy.pi_test(5)''',setup='import example_cy', number=100)
# cpi = timeit.timeit('''example_cy.cpi_test(5)''',setup='import example_cy', number=100)
# py = timeit.timeit('''example.test(5)''',setup='import example', number=100)


# print('Test is %f' % cy)
# print('Pi_Test is %f' % (pi))
# print('cPi_Test is %f' % (cpi))
# print('Py_Test is %f' % (py))

from ase import *
from cmath import pi

import sys
import numpy as np
sys.path.append('/home/sanantoniochili/Desktop/PhD/Scripts/Switch_Implementation/gradients_implementation/pysrc')
from example import Potential
from potential import Coulomb as PCoulomb
from example import Coulomb

atoms = Atoms("SrTiO3",

				  cell=[[4.00, 0.00, 0.00],
						[0.00, 4.00, 0.00],
						[0.00, 0.00, 4.00]],

				  positions=[[0, 0, 0],
							 [2, 2, 2],
							 [0, 2, 2],
							 [2, 0, 2],
							 [2, 2, 0]],
				  pbc=True)

charge_dict = {
	'O' : -2.,
	'Sr':  2.,
	'Ti':  4.,
	'Cl': -1.,
	'Na':  1.,
	'S' : -2.,
	'Zn':  2.
}

vects = np.array(atoms.get_cell())
volume = abs(np.linalg.det(vects))
N = len(atoms.get_positions())
accuracy = 0.00001  # Demanded accuracy of terms 
alpha = N**(1/6) * pi**(1/2) / volume**(1/3)
real_cut = (-np.log(accuracy))**(1/2)/alpha
recip_cut = 2*alpha*(-np.log(accuracy))**(1/2)

p = Coulomb()
pp = PCoulomb()

p.set_parameters(alpha=alpha, real_cut_off=real_cut,
						recip_cut_off=recip_cut, 
						chemical_symbols=np.array(atoms.get_chemical_symbols()),
						charge_dict=charge_dict, N=N,
						filename=None)
pp.set_parameters(alpha=alpha, real_cut_off=real_cut,
						recip_cut_off=recip_cut, 
						chemical_symbols=atoms.get_chemical_symbols(),
						charge_dict=charge_dict,
						filename=None)
pos = atoms.positions

import timeit
pt = timeit.timeit('''p.calc_real(pos,vects,N)''', globals=globals(), number=1)
ppt = timeit.timeit('''print(sum(sum(pp.calc_complete(N,pp.calc_real(pos,vects,N)))))''', globals=globals(), number=1)

print('Cython: %f' % pt)
print('Python: %f' % ppt)

# p = Coulomb()
# symbs = np.array(["Sr","Ti"], dtype=np.dtype("U2"))
# chs = {"Sr":1, "Ti":2}
# p.set_parameters(0,0,0,symbs,2,chs)