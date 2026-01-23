#!/usr/bin/env python3

import elphmod.models.graphene
import numpy as np
import subprocess

elphmod.misc.verbosity = 0

kT = 0.025 / elphmod.misc.Ry
nc = 12
n = nc * nc

el, ph, elph, elel = elphmod.models.graphene.create(rydberg=True)

e = np.linalg.eigvalsh(el.supercell(nc, nc).H())
mu = elphmod.occupations.find_Fermi_level(n, e, kT)

res = float(subprocess.check_output(['./elphy', str(nc)]))

ref = elphmod.diagrams.grand_potential(e - mu, kT) + mu * n

assert np.allclose(res, ref)
