# (c) 2017 Gregor Mitscha-Baude
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

params = dict(
    h = 0.5,
    Nmax = 3e4,
    bV = 0.,
)
name = "pot-ahem"
if not fields.exists(name, **params):
    setup = model.Setup(**params)
    pb, pnps = model.solve(setup)
    v = pnps.solutions()[0]
    fields.save_functions(name, params, v=v)

fun, mesh = fields.get_functions(name, **params)
f = fun["v"]

def F(x, z):
    if x>=0:
        return f([x, z])
    else:
        return f([-x, z])

R, Htop, Hbot = 7, 2, 12
N = 50
Nx, Ny = N+1, 2*N + 1

Y, X = np.mgrid[-Hbot:Htop:Ny*1j, -R:R:Nx*1j]
U = np.zeros((Ny,Nx))

for y in range(Ny):
    for x in range(Nx):
        U[y][x] = F(X[y][x], Y[y][x])

fig, ax = plt.subplots(figsize=(8, 6), num="pot")

tr, vertex_values = mesh2triang(mesh)
zz = vertex_values(f)
pc = ax.tripcolor(tr, zz, cmap=cm.coolwarm_r)
#pc = plt.pcolor(X, Y, U, cmap=cm.coolwarm_r) #, vmin=0, vmax=1)
plt.colorbar(pc)
#plot_polygon(ax, pugh.polygon(diamPore=6., rmem=13))
plt.xlim(-R, R)
plt.ylim(-Hbot, Htop)
minus = lambda l: [[-x[0], x[1]] for x in l]

p = get_pore()
p.protein.plot("-k")
Polygon(minus(p.protein.nodes)).plot("-k")

from folders import FIGDIR
from nanopores import savefigs
#savefigs("potential", DIR=FIGDIR + "/ahem")

#dolfin.plot(v, backend="matplotlib", cmap=cm.viridis)
plt.show()