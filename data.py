def put_model(filename, el, ph, elph, A, kT, n, nspin=2, eps=1e-10):
    Ry2Ha = 0.5

    A = elphmod.bravais.supercell(*A)[1]

    Ri = list(map(tuple, el.R))
    Rj = list(map(tuple, ph.R))
    Rk = list(map(tuple, elph.Rk))
    Rl = list(map(tuple, elph.Rk))

    R = sorted(set(Ri + Rj + Rk + Rl))

    hoppings = []

    for ri in range(len(Ri)):
        i = R.index(Ri[ri])
        for a in range(el.size):
            for b in range(el.size):
                t = el.data[ri, a, b].real * Ry2Ha
                if abs(t) > eps:
                    hoppings.append((i, a, b, t))

    springs = []

    for rj in range(len(Rj)):
        j = R.index(Rj[rj])
        for x in range(ph.size):
            for y in range(ph.size):
                k = ph.data[rj, x, y].real * Ry2Ha
                if abs(k) > eps:
                    springs.append((j, x, y, k))

    couplings = []

    for rk in range(len(Rk)):
        k = R.index(Rk[rk])
        for z in range(ph.size):
            for rl in range(len(Rl)):
                l = R.index(Rl[rl])
                for c in range(el.size):
                    for d in range(el.size):
                        g = elph.data[rk, z, rl, c, d].real * Ry2Ha
                        if abs(g) > eps:
                            couplings.append((k, z, l, c, d, g))

    with open(filename, 'w') as dat:
        dat.write(f'{kT * Ry2Ha}\n')
        dat.write(f'{n}\n')
        dat.write(f'{elph.el.size}\n')
        dat.write(f'{nspin}\n')

        for i in range(3):
            dat.write('%2d %2d %2d\n' % tuple(A[i]))

        for i in range(3):
            dat.write('%12.9f %12.9f %12.9f\n' % tuple(elph.ph.a[i]))

        dat.write(f'{elph.ph.nat}\n')

        for i in range(elph.ph.nat):
            dat.write('%2s %12.9f %12.9f %12.9f 0 0 0\n'
                % (elph.ph.atom_order[i], *elph.ph.r[i]))

        dat.write(f'{len(R)}\n')

        for r in R:
            dat.write('%2d %2d %2d\n' % r)

        dat.write(f'{len(hoppings)}\n')

        for t in sorted(hoppings):
            dat.write('%d %d %d %12.9f\n' % t)

        dat.write(f'{len(springs)}\n')

        for k in sorted(springs):
            dat.write('%d %d %d %12.9f\n' % k)

        dat.write(f'{len(couplings)}\n')

        for g in sorted(couplings):
            dat.write('%d %d %d %d %d %12.9f\n' % g)

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

    put_model('input.dat', el, ph, elph, A, kT, n)
    driver.save('driver.pickle')
