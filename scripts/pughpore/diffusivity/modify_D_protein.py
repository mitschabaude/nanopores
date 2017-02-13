from matplotlib import pyplot as plt
import os
from nanopores.tools import fields, smooth
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
#plt.show()
