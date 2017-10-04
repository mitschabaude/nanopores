from matplotlib import pyplot as plt
import numpy as np
import os
from nanopores.tools import fields, smooth, savefigs
from nanopores.models.diffusion_interpolation import Dn_plane, Dt_plane
fields.set_dir(os.path.expanduser("~") + "/Dropbox/nanopores/fields")

rMolecule = 2.0779
data = fields.get_fields("pugh_diff3D_cross", bulkbc=True, rMolecule=rMolecule)
R = 3.
x = [R - z[0] for z in data["x"]]
data, x = fields._sorted(data, x)

Dxx_ = [D[0][0] for D in data["D"]]
Dyy_ = [D[1][1] for D in data["D"]]
Dzz_ = [D[2][2] for D in data["D"]]

Dxx = smooth(smooth(Dxx_, 3), 5)
Dyy = smooth(smooth(Dyy_, 3), 5)
Dzz = smooth(smooth(Dzz_, 3), 5)

def matrix(d):
    return [[d[0], 0., 0.], [0., d[1], 0.], [0., 0., d[2]]]

data = dict(x=x, D=map(matrix, zip(Dxx, Dyy, Dzz)))
if not fields.exists("pugh_diff_pore", rMolecule=rMolecule):
    print "SAVING..."
    fields.save_fields("pugh_diff_pore", dict(rMolecule=rMolecule), **data)
    fields.update()

plt.figure()
plt.plot(x, Dxx_, "o:b")
plt.plot(x, Dxx, "-b", label=r"$D_x$")
plt.plot(x, Dyy_, "o:r")
plt.plot(x, Dyy, "-r", label=r"$D_y$")
plt.plot(x, Dzz_, "o:g")
plt.plot(x, Dzz, "-g", label=r"$D_z$")

plt.xlabel('x distance from pore wall [nm]')
plt.ylabel('diffusivity relative to bulk')
plt.legend(loc='lower right')
plt.tight_layout()

plt.figure()
xlin = np.linspace(rMolecule+1e-3, 3., 100)
dn = [Dn_plane(t, rMolecule, N=20) for t in xlin]
dx = [d*Dxx[-1]/dn[-1] for d in dn]

dt = [Dt_plane(t, rMolecule) for t in xlin]
dz = [d*Dzz[-1]/dt[-1] for d in dt]
dy = [d*Dyy[-1]/dt[-1] for d in dt]

plt.plot(xlin, dx, "-b")
plt.plot(xlin, dy, "-g")
plt.plot(xlin, dz, "-r")

plt.plot(x[::2], Dxx[::2], "ob", label=r"$D_{xx}$")
plt.plot(x[::2], Dyy[::2], "sg", label=r"$D_{yy}$")
plt.plot(x[::2], Dzz[::2], ".r", label=r"$D_{zz}$")

plt.axvline(x=rMolecule, linestyle="--", color="#666666")
plt.annotate("Protein radius", (rMolecule, 0.12),
                 xytext=(rMolecule + 0.1, 0.12-0.002), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))

plt.xlim(2., 3.)
plt.ylim(0., 0.13)
plt.gcf().set_size_inches(3.2, 3.2)
plt.xlabel("x distance from pore wall [nm]")
plt.ylabel("Rel. diffusivity")
plt.legend(loc='lower right')
plt.tight_layout()

HOME = os.path.expanduser("~")
savefigs("pugh_Dproteins", os.path.join(HOME, "Dropbox", "nanopores", "figures"))


