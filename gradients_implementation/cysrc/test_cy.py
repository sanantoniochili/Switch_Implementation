# cy = timeit.timeit('''example_cy.test(5)''',setup='import example_cy',number=100)
# pi = timeit.timeit('''example_cy.pi_test(5)''',setup='import example_cy', number=100)
# cpi = timeit.timeit('''example_cy.cpi_test(5)''',setup='import example_cy', number=100)
# py = timeit.timeit('''example.test(5)''',setup='import example', number=100)


# print('Test is %f' % cy)
# print('Pi_Test is %f' % (pi))
# print('cPi_Test is %f' % (cpi))
# print('Py_Test is %f' % (py))

import sys
import numpy as np
sys.path.append('/home/sanantoniochili/Desktop/PhD/Scripts/Switch_Implementation/gradients_implementation/pysrc')
from example import Potential
from potential import Potential as PPot

p = Potential()
pp = PPot()
cutoff = 1
vects = np.array([[1,1,1],[1,1,1],[1,1,1]], dtype=np.double)

import timeit
pt = timeit.timeit('''shifts = np.array(p.get_shifts(cutoff, vects));print(shifts)''', globals=globals(), number=1)
ppt = timeit.timeit('''shifts = np.array(pp.get_shifts(cutoff, vects));print(shifts)''', globals=globals(), number=1)

print('Cython Potential: %f' % pt)
print('Python Potential: %f' % ppt)