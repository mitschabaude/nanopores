# (c) 2017 Gregor Mitscha-Baude
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
from matplotlib import cm
import matplotlib.tri as tri
import nanopores.models.nanopore as model
import dolfin
from nanopores.tools import fields
from nanopores.geometries.alphahempoly import poly
from nanopores.geometries.alphahem import default, get_pore

from nanopores.tools.polygons import Polygon

def mesh2triang(mesh):
    x = mesh.coordinates()[:,0]
    y = mesh.coordinates()[:,1]
    triangles = mesh.cells()
    notx0 = x>0.
    x2 = -x[notx0]
    y2 = y[notx0]
    xx = np.concatenate([x, x2])
    yy = np.concatenate([y, y2])
    N0 = x.shape[0]
    triangles2 = np.array(triangles)
    # TODO: slow
    for i, j in enumerate(np.where(notx0)[0]):
        triangles2[triangles2==j] = i+N0
    tt = np.concatenate([triangles, triangles2])

    def vertex_values(u):
        z = u.compute_vertex_values(mesh)
        zz = np.zeros(xx.shape)
        zz[:N0] = z
        zz[N0:] = z[notx0]
        return zz

    return tri.Triangulation(xx, yy, tt), vertex_values

minus = lambda l: [[-x[0], x[1]] for x in l]

def plot_polygon(p, lw=1):
    p.plot("-k", lw=lw)
    Polygon(minus(p.nodes)).plot("-k", lw=lw)

params = dict(
    h = 0.5,
    Nmax = 3e4,
    bV = 0.0,
)
name = "pot-ahem"
if not fields.exists(name, **params):
    setup = model.Setup(**params)
    pb, pnps = model.solve(setup)
    v = pnps.solutions()[0]
    fields.save_functions(name, params, v=v)
    fields.update()

fun, mesh = fields.get_functions(name, **params)
f = fun["v"]

def F(x, z):
    if x>=0:
        return f([x, z])
    else:
        return f([-x, z])

R, Htop, Hbot = 7, 2, 12
fig = plt.figure(figsize=(2.6, 2.15))
ax = plt.axes(xlim=(-R, R), ylim=(-Hbot, Htop))

tr, vertex_values = mesh2triang(mesh)
zz = 1000*vertex_values(f)

minmax = max(abs(zz))
#plt.triplot(tr)
#pc = plt.tripcolor(tr, zz, cmap=cm.bwr)
pc = plt.tripcolor(tr, zz, cmap=cm.bwr, vmin=-minmax, vmax=minmax)
#pc = plt.pcolor(X, Y, U, cmap=cm.coolwarm_r) #, vmin=0, vmax=1)
cb = plt.colorbar(pc)
cb.ax.set_ylabel("Electric potential [mV]")
#plot_polygon(ax, pugh.polygon(diamPore=6., rmem=13))
#plt.xlim(-R, R)
#plt.ylim(-Hbot, Htop)

p = get_pore()
plot_polygon(p.protein, lw=.5)
plot_polygon(p.membrane, lw=.5)

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

from folders import FIGDIR_HOWORKA
from nanopores import savefigs
savefigs("potential", FIGDIR_HOWORKA + "/ahem", ending=".pdf")

#dolfin.plot(v, backend="matplotlib", cmap=cm.viridis)
plt.show()