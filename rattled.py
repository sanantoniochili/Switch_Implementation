from ase import *
rom ase.io import read, write
from ase.visualize import view
from ase.calculators.gulp import GULP


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


''' Perturb atoms '''
atoms.rattle(stdev=0.1)

 ase.io.write(filename, images, format=None, parallel=True, append=False, **kwargs)[source]