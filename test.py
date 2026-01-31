#!/usr/bin/env python3

from input import *

import numpy as np
import subprocess

driver = elphmod.md.Driver(elph, kT, 'fd', n, supercell=A, unscreen=False)

res = np.array([float(x.rstrip(';'))
    for x in subprocess.check_output(['./elphy', 'input.dat', 'input.xyz'],
        universal_newlines=True).split() if '.' in x])

driver.from_xyz('input.xyz')

energy = driver.free_energy(show=False)
forces = -driver.jacobian(show=False)

ref = np.insert(forces, 0, energy)

assert np.allclose(res, ref * Ry2Ha)
