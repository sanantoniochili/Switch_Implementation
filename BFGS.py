from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
import numpy as np
import tkinter

d = 0.9575
t = np.pi / 180 * 104.51
water = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)],
              calculator=EMT())
# dyn = BFGS(water, trajectory='H2O.traj')
dyn = BFGS(atoms=water, trajectory='qn.traj', restart='qn.pckl')
dyn.run(fmax=0.05)

