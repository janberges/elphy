import elphmod

driver = elphmod.md.Driver.load('driver.pickle')

with open('input.dat') as lines:
    for _ in range(4):
        next(lines)

    driver.elph.ph.r *= 1.0 + float(next(lines))

driver.plot(interactive=True, scale=10, pause=0.1)
driver.from_xyz('ipi.pos_0.xyz')
