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

    parameters = dict(kT=0.0019, n=2.0, supercell=(12, (6, 12, 0)))
    strain = 0.3

    elph.export(dat, strain=strain, **parameters)

    el.data *= 1 - elphmod.models.graphene.beta * strain
    ph.a *= 1 + strain
    ph.r *= 1 + strain

    driver = elphmod.md.Driver(elph, f='fd', unscreen=False, **parameters)

elif model == 'TaS2':
    import elphmod.models.tas2

    el, ph, elph = elphmod.models.tas2.create(rydberg=True, divide_mass=False)

    driver = elphmod.md.Driver(elph, kT=0.005, f='fd', n=1.0,
        nk=(12, 12), nq=(2, 2), supercell=(9, 9), kT0=0.02, f0='mv', export=dat)
else:
    elphmod.MPI.info(f'Usage: python3 {sys.argv[0]} '
        '[(graphene|TaS2) [<data file> [<xyz file>]]]', error=True)

def run(radius):
    out = subprocess.check_output(f'./elphy {dat} {xyz} {radius}'.split(),
        universal_newlines=True).split('\n')

    energy = float(out[1].split()[-1].split('=')[1].strip('"'))
    forces = np.array([float(x) for line in out[2:] for x in line.split()[-3:]])

    driver.from_xyz(xyz)

    energy_elphmod = 0.5 * driver.free_energy(show=False)
    forces_elphmod = -0.5 * driver.jacobian(show=False)

    return energy, forces, energy_elphmod, forces_elphmod

e0a, f0a, e0b, f0b = run(0.0)
e1a, f1a, e1b, f1b = run(0.1)

ok = np.allclose(f1a - f0a, f1b - f0b) and np.allclose(e1a - e0a, e1b - e0b)

elphmod.MPI.info(f'elphmod and elphy {"" if ok else "DO NOT "}agree!',
    error=not ok)
