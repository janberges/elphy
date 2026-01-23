#!/usr/bin/env python3

import elphmod.models.graphene

elphmod.misc.verbosity = 0

eps = 1e-10

el, ph, elph, elel = elphmod.models.graphene.create(rydberg=True,
    divide_mass=False)

def put_model(filename, model):
    elements = []

    for r in range(len(model.R)):
        for a in range(model.size):
            for b in range(model.size):
                c = model.data[r, a, b].real
                if abs(c) > eps:
                    elements.append((r, a, b, c))

    with open(filename, 'w') as dat:
        dat.write('%d %d %d\n' % (model.size, len(model.R), len(elements)))

        for r in model.R:
            dat.write('% d % d % d\n' % tuple(r))

        for t in sorted(elements):
            dat.write('%d %d %d % .9f\n' % t)

def put_coupling(filename, model):
    elements = []

    Rg = list(map(tuple, model.Rg))
    Rk = list(map(tuple, model.Rk))

    R = sorted(set(Rg + Rk))

    rg = [R.index(r) for r in Rg]
    rk = [R.index(r) for r in Rk]

    for g in range(len(model.Rg)):
        for x in range(model.ph.size):
            for k in range(len(model.Rk)):
                for a in range(model.el.size):
                    for b in range(model.el.size):
                        c = model.data[g, x, k, a, b].real
                        if abs(c) > eps:
                            elements.append((rg[g], x, rk[k], a, b, c))

    with open(filename, 'w') as dat:
        dat.write('%d %d %d %d\n' % (model.ph.size, model.el.size,
            len(R), len(elements)))

        for r in R:
            dat.write('% d % d % d\n' % tuple(r))

        for t in sorted(elements):
            dat.write('%d %d %d %d %d % .9f\n' % t)

put_model('el.dat', el)
put_model('ph.dat', ph)
put_coupling('elph.dat', elph)
