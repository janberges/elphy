#!/usr/bin/env python3

import elphmod.models.graphene
import numpy as np
import subprocess

elphmod.misc.verbosity = 0

nc = 12

el, ph, elph, elel = elphmod.models.graphene.create()

res = np.array(list(map(float, subprocess.check_output(['./elphy', str(nc)],
    universal_newlines=True).split())))

ref = np.linalg.eigvalsh(el.supercell(nc, nc).H())

assert np.allclose(res, ref)
