from ase import *
from ase.visualize import view
from ase.io import read as aread

import numpy as np
import pandas as pd
import sys,os
from pathlib import Path

DATAPATH = "../../Data/material/ICSD_SrTiO3"
# folder = input("Insert folder with cifs to use:\n")

for file in Path(DATAPATH).rglob('*.cif'):
	atoms = aread(file)
	symbols = atoms.get_chemical_symbols()
	# dists = atoms.get_all_distances()

	# print(file.name)
