import timeit

# cy = timeit.timeit('''example_cy.test(5)''',setup='import example_cy',number=100)
# pi = timeit.timeit('''example_cy.pi_test(5)''',setup='import example_cy', number=100)
# cpi = timeit.timeit('''example_cy.cpi_test(5)''',setup='import example_cy', number=100)
# py = timeit.timeit('''example.test(5)''',setup='import example', number=100)


# print('Test is %f' % cy)
# print('Pi_Test is %f' % (pi))
# print('cPi_Test is %f' % (cpi))
# print('Py_Test is %f' % (py))

pt = timeit.timeit('''p.set_parameters()''', setup='from example_cy import Potential;p = Potential()',number=1)
ppt = timeit.timeit('''pp.set_parameters()''', setup='from potential import Potential as PPot;pp=PPot()',number=1)

print('Cython Potential: %f' % pt)
print('Python Potential: %f' % ppt)