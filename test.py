#!/usr/bin/env python3

import elphmod.models.graphene
import numpy as np
import subprocess

elphmod.misc.verbosity = 0

nc = 12

elph = elphmod.models.graphene.create(rydberg=True, divide_mass=False)[2]
elph.data[...] = 0.0

driver = elphmod.md.Driver(elph, kT=0.0019, f='fd', n=0.5 * elph.el.size,
    supercell=(nc, nc), unscreen=False)

driver.random_displacements()

with open('u.dat', 'w') as dat:
    for u in driver.u:
        dat.write('% .9f\n' % u)

res = float(subprocess.check_output(['./elphy', str(nc)]))

ref = driver.free_energy(show=False)

assert np.allclose(res, ref)
