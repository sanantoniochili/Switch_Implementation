from ase import *
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
# view(atoms)

# set calculator
calc = GULP(keywords='conp', library='buck.lib')
atoms.calc = calc

# calculate energy
E= atoms.get_potential_energy()
print(E)

# perturb atoms
atoms.rattle(stdev=0.1)
# view(atoms)

# optimisation
calc.set(keywords='conp opti')
opt = calc.get_optimizer(atoms)
opt.run()
E = atoms.get_potential_energy()
print(E)