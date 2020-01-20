from sympy import *
import numpy as np

x, y = symbols('x y')
A, B, C = symbols('A B C')
q_1, q_2, p, e = symbols('q_1 q_2 p e')


f = A*(exp(-B*sqrt(x**2 + y**2))) - C/((x**2 + y**2)**3) + \
    (q_1*q_2/(4*p*e)) * (1/sqrt(x**2 + y**2))

deriv = f.diff(x)
print(deriv.evalf(subs={x: 1, y: 1}))