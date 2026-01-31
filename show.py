from data import *

elph.clear()
driver = elphmod.md.Driver(elph, kT, 'fd', n, supercell=A, unscreen=False)

driver.plot(interactive=True, scale=10, pause=0.1)
driver.from_xyz("ipi.pos_0.xyz")
