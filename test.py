import elphmod
import numpy as np
import subprocess

res = np.array([float(x.rstrip(';'))
    for x in subprocess.check_output(['./elphy', 'input.dat', 'input.xyz'],
        universal_newlines=True).split() if '.' in x])

driver = elphmod.md.Driver.load('driver.pickle')
driver.from_xyz('input.xyz')

energy = driver.free_energy(show=False)
forces = driver.F0 - driver.jacobian(show=False)

ref = 0.5 * np.insert(forces, 0, energy)

assert np.allclose(res, ref)
