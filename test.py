#!/usr/bin/env python3

import elphmod.models.graphene
import numpy as np
import subprocess

elphmod.misc.verbosity = 0

nc = 12
n = nc * nc

el, ph, elph, elel = elphmod.models.graphene.create()

e = np.linalg.eigvalsh(el.supercell(nc, nc).H())
mu = elphmod.occupations.find_Fermi_level(n, e)

res = float(subprocess.check_output(['./elphy', str(nc)]))

ref = elphmod.diagrams.grand_potential(e - mu) + mu * n

assert np.allclose(res, ref)
