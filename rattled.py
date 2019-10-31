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

# perturb atoms
atoms.rattle(stdev=0.1)

''' Set calculator. unfix keyword is important; 
by default GULP will set the first derivative 
of the first atom in the list to zero,
since this atom is not usually allowed to move 
during relaxations in GULP '''
calc = GULP(keywords='opti conj conp full nosymm',
            library='buck.lib')
            # trajectory=traj)
atoms.set_calculator(calc)


''' Calculate energy '''
E = atoms.get_potential_energy()
print(E)
