import elphmod
import numpy as np
import subprocess
import sys

indat = sys.argv[1] if len(sys.argv) > 1 else 'input.dat'
model = sys.argv[2] if len(sys.argv) > 2 else 'graphene'
units = sys.argv[3] if len(sys.argv) > 3 else 'Ha'

def error():
    elphmod.MPI.info(f'Usage: python3 {sys.argv[0]} '
        '[<data file> [(graphene|TaS2|chain|Be) [(Ha|Ry|eV)]]]', error=True)

if units == 'Ha':
    econv = 0.5
    lconv = 1.0
    mconv = 2.0
elif units == 'Ry':
    econv = 1.0
    lconv = 1.0
    mconv = 1.0
elif units == 'eV':
    econv = elphmod.misc.Ry
    lconv = elphmod.misc.a0
    mconv = 2 * elphmod.misc.meSI * (1e-10 / 1e-15) ** 2 / elphmod.misc.eVSI
    # chosen such that time is measured in femtoseconds
else:
    error()

create = dict(rydberg=True, divide_mass=False)
export = dict(export=indat, econv=econv, lconv=lconv, mconv=mconv)

if model == 'graphene':
    import elphmod.models.graphene

    el, ph, elph, elel = elphmod.models.graphene.create(**create)

    parameters = dict(kT=0.0019, n=2.0, supercell=(12, (6, 12, 0)))
    strain = 0.3

    export['filename'] = export.pop('export')
    elph.export(strain=strain, **parameters, **export)

    el.data *= 1 - elphmod.models.graphene.beta * strain
    ph.a *= 1 + strain
    ph.r *= 1 + strain

    driver = elphmod.md.Driver(elph, f='fd', unscreen=False, **parameters)

elif model == 'TaS2':
    import elphmod.models.tas2

    el, ph, elph = elphmod.models.tas2.create(**create)

    driver = elphmod.md.Driver(elph, kT=0.005, f='fd', n=1.0, supercell=(9, 9),
        nk=(12, 12), nq=(2, 2), kT0=0.02, f0='mv', **export)

elif model == 'chain':
    import elphmod.models.chain

    el, ph, elph = elphmod.models.chain.create(**create)

    driver = elphmod.md.Driver(elph, kT=1e-3, f='fd', n=1.0, supercell=(23,),
        unscreen=False, **export)

elif model == 'Be':
    import elphmod.models.be

    el, ph, elph = elphmod.models.be.create(**create)

    driver = elphmod.md.Driver(elph, kT=1e-3, f='fd', n=2.0, supercell=(12, 12),
        unscreen=False, **export)
else:
    error()

def run(radius):
    out = subprocess.check_output(f'./elphy {indat} 1 {radius}'.split(),
        universal_newlines=True).split('\n')

    energy = float(out[1].split()[-1].split('=')[1].strip('"'))
    rows = np.array([list(map(float, line.split()[1:])) for line in out[2:-1]])

    driver.u = (rows[:, :3] / lconv - driver.elph.ph.r).ravel()
    forces = rows[:, 3:].ravel()

    energy_elphmod = econv * driver.free_energy(show=False)
    forces_elphmod = -econv / lconv * driver.jacobian(show=False)

    return energy, forces, energy_elphmod, forces_elphmod

e0a, f0a, e0b, f0b = run(0.0)
e1a, f1a, e1b, f1b = run(0.1)

ok = np.allclose(f0a, f0b) and np.allclose(f1a, f1b)

if model == 'graphene':
    ok = ok and np.allclose(e1a - e0a, e1b - e0b)
else:
    ok = ok and np.allclose(e0a, e0b) and np.allclose(e1a, e1b)

elphmod.MPI.info(f'elphmod and elphy {"" if ok else "DO NOT "}agree!',
    error=not ok)

tconv = np.sqrt(mconv / econv) * lconv

if elphmod.MPI.comm.rank == 0:
    with open(indat, 'a') as data:
        data.write(f"""Snippets for MD makefile target:
./elphy {indat} -1 {0.04 * lconv}
./elphy {indat} 1001:10 {20 * tconv} {0.0004 / tconv}:0 0
""")
