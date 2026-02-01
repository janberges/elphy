import elphmod.models.graphene

elphmod.misc.verbosity = 0

socket = 'localhost:31415'
Ry2Ha = 0.5
kT = 0.0019
umax = 0.1
n = 2.0
nspin = 2
A = elphmod.bravais.supercell(12, (6, 12, 0))[1]
eps = 1e-10

el, ph, elph, elel = elphmod.models.graphene.create(rydberg=True,
    divide_mass=False)

elph.data *= 1.5 # otherwise the system is stable

def put_model(filename):
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
        dat.write(f'{socket}\n')

        dat.write(f'{kT * Ry2Ha}\n')
        dat.write(f'{umax}\n')
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
    put_model('input.dat')
