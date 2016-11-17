# (c) 2016 Gregor Mitscha-Baude
import diffusion
import nanopores.tools.solvers as solvers
import nanopores.models.pughpore as pugh
import numpy as np
from matplotlib import pyplot as plt
from itertools import product
import dolfin
import nanopores
#import nanopores.tools.fields as fields
from folders import fields

data = fields.get_fields("pugh_diff3D_test")

def _sorted(data, key):
    I = sorted(range(len(key)), key=lambda k: key[k])
    return {k: [data[k][i] for i in I] for k in data}, [key[i] for i in I]

x = [z[0] for z in data["x"]]
data, x = _sorted(data, x)
dstr = ["x", "y", "z"]
for i in range(3):
    Dxx = [D[i][i] for D in data["D"]]
    style = "s-"
    plt.plot(x, Dxx, style, label=r"$D_{%s%s}$" % (dstr[i], dstr[i]))
plt.legend(loc="best")
plt.show()