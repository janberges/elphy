#!/usr/bin/env python3

import elphmod.models.graphene
import numpy as np
import subprocess

elphmod.misc.verbosity = 0

nc = 12

elph = elphmod.models.graphene.create(rydberg=True, divide_mass=False)[2]

driver = elphmod.md.Driver(elph, kT=0.0019, f='fd', n=0.5 * elph.el.size,
    supercell=(nc, nc), unscreen=False)

driver.random_displacements()

with open('u.dat', 'w') as dat:
    for u in driver.u:
        dat.write('% .9f\n' % u)

res = np.array(list(map(float, subprocess.check_output(['./elphy', str(nc)],
    universal_newlines=True).split())))

energy = driver.free_energy(show=False)
forces = driver.jacobian(show=False)

ref = np.insert(forces, 0, energy)

assert np.allclose(res, ref)
