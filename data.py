import copy
import elphmod
import numpy as np

def put_model(filename, elph, A, kT, n, nspin=2, strain=0.0, eps=1e-10):
    """Export electron-phonon model from `elphmod` to input file for `elphy`.

    Parameters
    ----------
    filename : str
        Name of input file for `elphy`.
    elph : obj
        Localized model for electron-phonon coupling.
    A : ndarray
        Supercell lattice vectors in units of primitive lattice vectors.
    kT : float
        Smearing temperature in Ry.
    n : float
        Number of electrons per primitive cell.
    nspin : int, default 2
        Number of spins per orbital.
    strain : int, default 0.0
        Isotropic strain.
    eps : float, defailt 1e-10
        Matrix-element threshold in Hartree atomic units.
    """
    Ry2Ha = 0.5

    A = elphmod.bravais.supercell(*A)[1]

    elph = copy.deepcopy(elph)

    elph.el.standardize(eps=eps / Ry2Ha / abs(elph.el.data).max())
    elph.ph.standardize(eps=eps / Ry2Ha / abs(elph.ph.data).max())
    elph.standardize(eps=eps / Ry2Ha / abs(elph.data).max())

    Ri = list(map(tuple, elph.el.R))
    Rj = list(map(tuple, elph.ph.R))
    Rk = list(map(tuple, elph.Rk))
    Rl = list(map(tuple, elph.Rk))

    R = sorted(set(Ri + Rj + Rk + Rl))

    hoppings = []

    for ri in range(len(Ri)):
        i = R.index(Ri[ri])
        for a in range(elph.el.size):
            for b in range(elph.el.size):
                t = elph.el.data[ri, a, b].real * Ry2Ha
                if t:
                    hoppings.append((i, a, b, t))

    springs = []

    for rj in range(len(Rj)):
        j = R.index(Rj[rj])
        for x in range(elph.ph.size):
            for y in range(elph.ph.size):
                k = elph.ph.data[rj, x, y].real * Ry2Ha
                if k:
                    springs.append((j, x, y, k))

    couplings = []

    for rk in range(len(Rk)):
        k = R.index(Rk[rk])
        for z in range(elph.ph.size):
            for rl in range(len(Rl)):
                l = R.index(Rl[rl])
                for c in range(elph.el.size):
                    for d in range(elph.el.size):
                        g = elph.data[rk, z, rl, c, d].real * Ry2Ha
                        if g:
                            couplings.append((k, z, l, c, d, g))

    with open(filename, 'w') as dat:
        dat.write(f'{kT * Ry2Ha}\n')
        dat.write(f'{n}\n')
        dat.write(f'{elph.el.size}\n')
        dat.write(f'{nspin}\n')
        dat.write(f'{strain}\n')

        for i in range(3):
            dat.write('%2d %2d %2d\n' % tuple(A[i]))

        for i in range(3):
            dat.write('%15.9f %15.9f %15.9f\n' % tuple(elph.ph.a[i]))

        dat.write(f'{elph.ph.nat}\n')

        for i in range(elph.ph.nat):
            dat.write('%2s %15.9f %15.9f %15.9f 0 0 0\n'
                % (elph.ph.atom_order[i], *elph.ph.r[i]))

        dat.write(f'{len(R)}\n')

        for r in R:
            dat.write('%2d %2d %2d\n' % r)

        f = f'% .{max(0, -int(np.floor(np.log10(2 * eps))))}f'

        dat.write(f'{len(hoppings)}\n')

        for t in sorted(hoppings):
            dat.write(f'%d %d %d {f}\n' % t)

        dat.write(f'{len(springs)}\n')

        for k in sorted(springs):
            dat.write(f'%d %d %d {f}\n' % k)

        dat.write(f'{len(couplings)}\n')

        for g in sorted(couplings):
            dat.write(f'%d %d %d %d %d {f}\n' % g)

if __name__ == '__main__':
    import sys

    model = sys.argv[1] if len(sys.argv) > 1 else 'none'

    if model == 'graphene':
        import elphmod.models.graphene

        A = 12, (6, 12, 0)
        kT = 0.0019
        n = 2.0

        el, ph, elph, elel = elphmod.models.graphene.create(rydberg=True,
            divide_mass=False)

        elph.data *= 1.5 # otherwise the system is stable

    elif model == 'TaS2':
        import elphmod.models.tas2

        A = 9, 9
        kT = 0.005
        n = 1.0

        el, ph, elph = elphmod.models.tas2.create(rydberg=True,
            divide_mass=False)

        driver = elphmod.md.Driver(elph, 0.02, 'mv', n, nk=(12, 12), nq=(2, 2))

        elphmod.ph.q2r(ph, nq=driver.nq, D_full=driver.C0, divide_mass=False)

    else:
        sys.exit(f'Usage: python3 {sys.argv[0]} graphene|TaS2')

    driver = elphmod.md.Driver(elph, kT, 'fd', n, supercell=A, unscreen=False)

    put_model('input.dat', elph, A, kT, n)
    driver.save('driver.pickle')
