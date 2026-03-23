import elphmod
import subprocess

subprocess.check_output('./elphy input.dat symmetric.xyz 0'.split())

typ, tau0 = next(elphmod.misc.read_xyz('symmetric.xyz'))

plot = elphmod.plot.AtomsPlot(tau0, typ)
plot.plot(interactive=True, scale=10, pause=0.1)

for typ, tau in elphmod.misc.read_xyz('ipi.pos_0.xyz'):
    plot.set_positions(tau)
    plot.set_displacements(tau - tau0)
    plot.update()
