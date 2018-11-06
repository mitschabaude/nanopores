# (c) 2018 Gregor Mitscha-Baude
from matplotlib import rcParams, rc
rcParams.update({
    "font.size" : 7,
    "axes.titlesize" : 7,
    "font.family" : "sans-serif",
    "font.sans-serif" : ["CMU Sans Serif"],
    "lines.linewidth" : 1,
    "lines.markersize" : 5,
})

import matplotlib.pyplot as plt
import numpy as np
import nanopores
import nanopores.geometries.pughpore as pughpore
from nanopores.models.pughpore import polygon as pughpolygon
from nanopores.models.pughpoints import plot_polygon
from nanopores.tools import fields
fields.set_dir_mega()
from nanopores.models.diffusion_interpolation import get_diffusivity

#    from nanopores.tools.utilities import uCross, RectangleMesh
#    from math import pi, sqrt
#
#    dparams = {2: dict(diamPore=6., diamDNA=2.5, Nmax=1.2e5, dim=2, r=0.11, h=.75,
#                       cheapest=False, Membraneqs=-.5),
#               3: dict(diamPore=6., Nmax=1e6, dim=3, r=0.11, h=2.0, cheapest=False)}

# obtain diffusivity field and project to x-z plane
params = dict(geoname = "pughcyl", dim=2, r=0.11, h=.5, Nmax=1e5,
              cheapest=False, Membraneqs=-0.2)
functions, mesh = get_diffusivity(**params)
#functions = get_pugh_diffusivity(**dparams[2])
#setup = pugh.Setup(dim=2, h=1., Nmax=1e5, x0=None, diffusivity="Dpugh2")
#setup.prerefine()
#pugh.set_D(setup)
#D3D = setup.phys.Dp[1, 1]
#print D3D([0.,0.])
D3D = functions["D"][0]

D0 = nanopores.D
def F(x, z):
    if x>=0:
        return D3D([x, z])/D0
    else:
        return D3D([-x, z])/D0
#D = uCross(u=D3D, axis=1, degree=1, dim=2)

# obtain 2D mesh where we will evaluate field
rx, ry = pughpore.params["R"], 0.5*pughpore.params["H"]
rx, ry = 15, 28
Nx, Ny = 201, 401
#mesh2D = RectangleMesh([-R,-H/2.], [R, H/2.], int(4*R), int(2*H))

Y, X = np.mgrid[-ry:ry:Ny*1j, -rx:rx:Nx*1j]
U = np.zeros((Ny,Nx))

for y in range(Ny):
    for x in range(Nx):
        U[y][x] = F(X[y][x], Y[y][x])

fig, ax = plt.subplots(figsize=(1.73, 1.9)) #, dpi=300)
pc = plt.pcolor(X, Y, U, cmap=plt.get_cmap("bone"), vmin=0, vmax=1)
plot_polygon(ax, pughpolygon(diamPore=6., rmem=15))
plt.xlim(-15, 15)
plt.ylim(-25, 28)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
bbox = ax.get_position()
l, b, h = bbox.x0, bbox.y0, bbox.height
w = 0.05
cbaxes = fig.add_axes([l - w - 0.07, b, w, h])
cb = plt.colorbar(pc, cax=cbaxes, ax=ax)    
cbaxes.set_ylabel("Rel. diffusivity") # (r"$D_{zz} / D_0$")
cb.set_ticks([0., 1.])
cb.set_ticklabels([0, 1])
cbaxes.yaxis.set_ticks_position('left')
cbaxes.yaxis.set_label_position('left')
cbaxes.yaxis.labelpad = -3.

import os
HOME = os.path.expanduser("~")
FIGDIR = os.path.join(HOME, "Dropbox", "Paper Howorka", "figures")
nanopores.savefigs("pugh/Dfield", FIGDIR)
plt.show()
