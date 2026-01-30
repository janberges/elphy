#!/usr/bin/env python3

import elphmod.models.graphene
import numpy as np
import subprocess

Ry2Ha = 0.5

elphmod.misc.verbosity = 0

elph = elphmod.models.graphene.create(rydberg=True, divide_mass=False)[2]

driver = elphmod.md.Driver(elph, kT=0.0019, f='fd', n=0.5 * elph.el.size,
    supercell=(12, 12), unscreen=False)

driver.random_displacements()
driver.to_xyz('input.xyz')

res = np.array([float(x.rstrip(';'))
    for x in subprocess.check_output(['./elphy', 'input.dat', 'input.xyz'],
        universal_newlines=True).split() if '.' in x])

energy = driver.free_energy(show=False)
forces = -driver.jacobian(show=False)

ref = np.insert(forces, 0, energy)

assert np.allclose(res, ref * Ry2Ha)
