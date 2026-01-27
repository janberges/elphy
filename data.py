#!/usr/bin/env python3

import elphmod.models.graphene

elphmod.misc.verbosity = 0

kT = 0.0019
n = 1.0
A = elphmod.bravais.supercell(12, 12)[1]
eps = 1e-10

el, ph, elph, elel = elphmod.models.graphene.create(rydberg=True,
    divide_mass=False)

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
            t = el.data[ri, a, b].real
            if abs(t) > eps:
                hoppings.append((i, a, b, t))

springs = []

for rj in range(len(Rj)):
    j = R.index(Rj[rj])
    for x in range(ph.size):
        for y in range(ph.size):
            k = ph.data[rj, x, y].real
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
                    g = elph.data[rk, z, rl, c, d].real
                    if abs(g) > eps:
                        couplings.append((k, z, l, c, d, g))

with open('model.dat', 'w') as dat:
    dat.write(f'{kT}\n')
    dat.write(f'{n}\n')
    dat.write(f'{elph.el.size}\n')

    for i in range(3):
        dat.write('%2d %2d %2d\n' % tuple(A[i]))

    dat.write(f'{elph.ph.size}\n')

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
