import elphmod
import numpy as np
import subprocess
import sys

model = sys.argv[1] if len(sys.argv) > 1 else 'graphene'
dat = sys.argv[2] if len(sys.argv) > 2 else 'input.dat'
xyz = sys.argv[3] if len(sys.argv) > 3 else 'input.xyz'

if model == 'graphene':
    import elphmod.models.graphene

    el, ph, elph, elel = elphmod.models.graphene.create(rydberg=True,
        divide_mass=False)

    elph.data *= 1.5 # otherwise the system is stable

    driver = elphmod.md.Driver(elph, kT=0.0019, f='fd', n=2.0,
        supercell=(12, (6, 12, 0)), unscreen=False, export=dat)

elif model == 'TaS2':
    import elphmod.models.tas2

    el, ph, elph = elphmod.models.tas2.create(rydberg=True, divide_mass=False)

    driver = elphmod.md.Driver(elph, kT=0.005, f='fd', n=1.0,
        nk=(12, 12), nq=(2, 2), supercell=(9, 9), kT0=0.02, f0='mv', export=dat)
else:
    elphmod.MPI.info(f'Usage: python3 {sys.argv[0]} '
        '[(graphene|TaS2) [<data file> [<xyz file>]]]', error=True)

res = np.array([float(x.rstrip(';'))
    for x in subprocess.check_output(f'./elphy {dat} {xyz} 0.1'.split(),
        universal_newlines=True).split() if '.' in x])

driver.from_xyz(xyz)

energy = driver.free_energy(show=False)
forces = driver.F0 - driver.jacobian(show=False)

ref = 0.5 * np.insert(forces, 0, energy)

assert np.allclose(res, ref)
