import elphmod
import sys

ref = sys.argv[1] if len(sys.argv) > 1 else 'symmetric.xyz'
xyz = sys.argv[2] if len(sys.argv) > 2 else 'ipi.pos_0.xyz'

typ, tau0 = next(elphmod.misc.read_xyz(ref))

plot = elphmod.plot.AtomsPlot(tau0, typ)
plot.plot(interactive=True, scale=10, pause=0.1)

for typ, tau in elphmod.misc.read_xyz(xyz):
    plot.set_positions(tau)
    plot.set_displacements(tau - tau0)
    plot.update()
