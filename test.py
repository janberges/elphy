import elphmod
import numpy as np
import subprocess
import sys

model = sys.argv[1] if len(sys.argv) > 1 else 'none'

if model == 'graphene':
    import elphmod.models.graphene

    el, ph, elph, elel = elphmod.models.graphene.create(rydberg=True,
        divide_mass=False)

    elph.data *= 1.5 # otherwise the system is stable

    driver = elphmod.md.Driver(elph, kT=0.0019, f='fd', n=2.0,
        supercell=(12, (6, 12, 0)), unscreen=False, export='input.dat')

elif model == 'TaS2':
    import elphmod.models.tas2

    el, ph, elph = elphmod.models.tas2.create(rydberg=True, divide_mass=False)

    driver = elphmod.md.Driver(elph, kT=0.005, f='fd', n=1.0, nk=(12, 12),
        nq=(2, 2), supercell=(9, 9), kT0=0.02, f0='mv', export='input.dat')
else:
    elphmod.MPI.info(f'Usage: python3 {sys.argv[0]} graphene|TaS2', error=True)

driver.save('driver.pickle')

res = np.array([float(x.rstrip(';'))
    for x in subprocess.check_output('./elphy input.dat input.xyz 0.1'.split(),
        universal_newlines=True).split() if '.' in x])

driver.from_xyz('input.xyz')

energy = driver.free_energy(show=False)
forces = driver.F0 - driver.jacobian(show=False)

ref = 0.5 * np.insert(forces, 0, energy)

assert np.allclose(res, ref)
