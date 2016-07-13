import dolfin
import nanopores
from math import sqrt
import colormaps as cm
import matplotlib
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import matplotlib.patches as patches
import numpy as np
from forcefield import geo, phys, Fel, Fdrag, params

#geo_name = "H_cyl_geo"
#geo_params = dict(
#x0 = None,
#rMolecule = 0.5,
#lcCenter = 0.05,
#lcMolecule = 0.05,
#)
#h = 8.
#nanopores.generate_mesh(h, geo_name, **geo_params)
#geo3D = nanopores.geo_from_name(geo_name, **geo_params)
#mesh = geo3D.mesh

mesh2D = nanopores.RectangleMesh([-10.,-10.], [10.,10.], 100, 100)
dolfin.plot(mesh2D)

dolfin.plot(Fel)
dolfin.plot(Fdrag)

Fel2, Fdrag2 = nanopores.convert2D(mesh2D, Fel, Fdrag)
dolfin.plot(Fel2)
dolfin.plot(Fdrag2)

Ny = 40
Nx = 50

bar, fig = plt.subplots(figsize=(5,4))
ax=plt.axes()

Y, X = np.mgrid[-10:10:Ny*1j, -10:10:Nx*1j] #-5:30 and -30:30
U = np.zeros((Ny,Nx))
V = np.zeros((Ny,Nx))
for y in range(Ny):
    for x in range(Nx):
        F=Fdrag2(X[y][x],Y[y][x])
        U[y][x] = F[0]
        V[y][x] = F[1]

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if a=="1.0":
        return r"$10^{{{}}}$".format(b)
    elif a=="0.0":
        return r"$0$"
    else:
        return r'${}\cdot10^{{{}}}$'.format(a,b)
formt = matplotlib.ticker.FuncFormatter(fmt)

strength = np.sqrt(U*U+V*V)
norm = matplotlib.colors.SymLogNorm(linthresh=1e-16, linscale=1.0, vmin=0., vmax=np.max(strength))
strm = plt.streamplot(X,Y,U,V,arrowsize=1.5, linewidth=1.5, density=1.5, cmap=cm.viridis, color=strength, norm=norm)
plt.colorbar(strm.lines, ticks=[0, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11], format=formt)
plt.xlabel('x [nm]') #, fontsize=20)
plt.ylabel('z [nm]') #, fontsize=20)

# plot pore polygon on top
X_How_2d = np.array([[1.,4.5],[2.5,4.5],[2.5,1.1],[10.,1.1],[10.,-1.1],[2.5,-1.1],[2.5,-4.5],[1.,-4.5]])
X_How_2dm = np.column_stack([-X_How_2d[:,0], X_How_2d[:,1]])

patch = patches.Polygon(X_How_2d, closed=True, facecolor="#eeeeee", linewidth=3., edgecolor="black")
patchm = patches.Polygon(X_How_2dm, closed=True, facecolor="#eeeeee", linewidth=3., edgecolor="black")
patch.set_zorder(10)
patchm.set_zorder(10)
ax.add_patch(patch)
ax.add_patch(patchm)

plt.show()

#nanopores.plot_cross_vector(Fel2, mesh2D, title="Fel")
#nanopores.plot_cross_vector(Fdrag2, mesh2D, title="Fdrag")
#dolfin.interactive()
#
#
#X_How_2d = np.array([[1.,4.5],[2.5,4.5],[2.5,1.1],[10.,1.1],[10.,-1.1],[2.5,-1.1],[2.5,-4.5],[1.,-4.5]])
#
#bar, fig = plt.subplots(figsize=(12,8))
#ax = plt.axes()
#
#def radius(*X):
#    return sqrt(sum(x**2 for x in X))
#
#leftend=15.
##x_mem=np.linspace(X_How_2d[18][0],leftend,100)
##y_mem=np.zeros(x_mem.shape[0])+X_How_2d[18][1]
##x_mem_2=-x_mem
#size=X_How_2d.shape[0]
#X=np.zeros(size+1)
#Y=np.zeros(size+1)
#for index in range(size):
#	X[index]=X_How_2d[index][0]
#	Y[index]=X_How_2d[index][1]
#X[size]=X[0]
#Y[size]=Y[0]
#X_2=-X
#
## whole domain: fac=0.1,p2=[
#
#axes=plt.gca()
#axes.set_ylim([-5,10])
#axes.set_xlim([-10,10])

