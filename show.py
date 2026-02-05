import elphmod

driver = elphmod.md.Driver.load('driver.pickle')
driver.plot(interactive=True, scale=10, pause=0.1)
driver.from_xyz('ipi.pos_0.xyz')
