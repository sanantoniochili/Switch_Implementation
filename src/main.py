# import torch
import sys
import argparse
import fileinput
from ase import *
from ase.io import read as aread

# x = torch.ones(2, 2, requires_grad=True)

# f = A*(exp(-B*sqrt(x**2 + y**2))) - C/((x**2 + y**2)**3) + \
#     (q_1*q_2/(4*p*e)) * (1/sqrt(x**2 + y**2))
# print(y)

# # vector-Jacobian product
# v = torch.tensor([[0.1, 1.0],[1, 0.0001]], dtype=torch.float)
# y.backward(v)

Atoms = aread("/users/phd/tonyts/Desktop/Data/RandomStart_Sr3Ti3O9/1.cif")
print(Atoms)