"tools to save/load/interpolate/modify forcefields given by node-value data"

"""
+ RBF interpolation
+ maybecalculate
+ vergleich 3D vs 2D
+ finite size vs. finite size interpolated
"""
import numpy as np
import math
import dolfin
from itertools import product
import nanopores
from nanopores.tools import fields
from nanopores.models import Howorka
#from plot_forcefield import porestreamlines

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.delaunay as dln
import matplotlib.ticker
import matplotlib.patches as patches
import colormaps as cm

nanopores.add_params(
    eps = 0.1,
)

def porestreamlines(polygon=None, rx=10., ry=10., Nx=100, Ny=100,
                    maxvalue=None, **fields):
    "streamlines plot of vector field around nanopore"  
    
    # interpolate on regular mesh symmetric w.r.t. center axis
    #mesh2D = nanopores.RectangleMesh([-rx-0.1,-ry-0.1], [rx+0.1,ry+0.1], Nx, Ny)
    #fields2 = nanopores.convert2D(mesh2D, *(fields.values()))    
    
    # prepare polygon and copy to left half
    settings = dict(closed=True, facecolor="#eeeeee", linewidth=3.,
                    edgecolor="black")
    if polygon is not None:
        polygon = np.array(polygon)
        polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])

    # prepare plots
    Ny += 1
    Nx += 1
    Y, X = np.mgrid[-ry:ry:Ny*1j, -rx:rx:Nx*1j]
    U = np.zeros((Ny,Nx))
    V = np.zeros((Ny,Nx))
    formt = matplotlib.ticker.FuncFormatter(fmt)
    ticks = [0] + [10**n for n in range(-16, -9)]
    
    # determine uniform color range from fields (maybe round to nearest 10-power)
    if maxvalue is None:
        maxvalue = max(dolfin.norm(F.vector(), "linf") for F in fields.values())
        #maxvalue = 10**int(np.log10(maxvalue))
    
    for i, F in enumerate(fields.values()):
        Fstr = fields.keys()[i]
        fig, ax = plt.subplots(figsize=(rx+1., ry), num=Fstr)        
        
        # fill array with function values
        for y in range(Ny):
            for x in range(Nx):
                f = F(X[y][x], Y[y][x])
                U[y][x] = f[0]
                V[y][x] = f[1]

        # streamplot with logarithmic scale
        strength = np.sqrt(U*U+V*V)
        norm = matplotlib.colors.SymLogNorm(linthresh=ticks[1], linscale=1.0,
                                            vmin=0., vmax=maxvalue)
        strm = plt.streamplot(X, Y, U, V, arrowsize=1., linewidth=1., density=3.,
                              cmap=cm.viridis, color=strength, norm=norm)
        plt.colorbar(strm.lines, ticks=ticks, format=formt)
        plt.xlabel('x [nm]') #, fontsize=20)
        plt.ylabel('z [nm]') #, fontsize=20)

        # plot pore polygon on top        
        if polygon is not None:
            patch = patches.Polygon(polygon, **settings)
            patchm = patches.Polygon(polygon_m, **settings)
            patch.set_zorder(10)
            patchm.set_zorder(10)
            ax.add_patch(patch)
            ax.add_patch(patchm)
            
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if a=="1.0":
        return r"$10^{{{}}}$".format(b)
    elif a=="0.0":
        return r"$0$"
    else:
        return r'${}\cdot10^{{{}}}$'.format(a,b)

# RBF, vectorized
def Gaussian(x, eps=1.):
    # x ... array of 3-vectors; x.shape = (..,3)
    d = len(x.shape)-1
    return np.exp(-np.sum(x**2, axis=d)/eps)    
def diff(x, y):
    # (xi), (yi) --> (xi-yj)i,j where xi, yj \in R^3
    # i.e. output has shape (N,M,3)
    return np.array(x)[:,None,:] - np.array(y)[None,:,:]

# build point-evaluation interpoland 
def interpoland(x, y, eps=1.):
    G = Gaussian(diff(x, x), eps)
    c = np.linalg.solve(G, y)
    def fv(x1): # vectorized interpoland
        return np.dot(Gaussian(diff(x1, x), eps), c)
    def f(x1): # single-eval interpoland
        return np.dot(Gaussian(x1 - x, eps), c)
    return f
    
# build discrete vector S1 interpoland on sample mesh
def S1(f, mesh, dim=1):
    V = dolfin.FunctionSpace(mesh, "CG", 1)
    if dim>1:
        V = dolfin.MixedFunctionSpace([V]*dim)
    f1 = dolfin.Function(V)
    value_shape = () if dim==1 else (dim,)
    class expr(dolfin.Expression):
        def eval(self, value, x):
            value[:] = f(x)
        def value_shape(self):
            return value_shape
    f1.interpolate(expr())
    return f1
    
# load force field and points
data = fields.get_fields("force3D", Qmol=-1, rMolecule=0.5)
X = data["x"]
F = data["F"]
xparams = fields.load_file("xforce")["params"]
Rx, Ry = xparams["Rx"], xparams["Ry"]

forces, mesh, params = nanopores.load_vector_functions("howorka2D_implicit")
Fimpl = forces["F"]
Fexp = [[1e-12*f[0], 1e-12*f[2]] for f in F]

# add additional points where implicit model will be used
def uniform_grid(xi, yi, h):
    Lx = xi[1]-xi[0]
    Ly = yi[1]-yi[0]
    m = int(math.ceil(Lx/h))
    n = int(math.ceil(Ly/h))
    x = np.linspace(xi[0], xi[1], m)
    y = np.linspace(yi[0], yi[1], n)
    X = list(product(x, y))
    return X
    
def scatter(X, **params):
    x = np.array([t[0] for t in X])
    y = np.array([t[1] for t in X])
    plt.scatter(x, y, **params)
    
from nanopores.geometries.H_cyl_geo.params_geo import r0, rMolecule, l0, r1
r = rMolecule
Ry0 = 7.
Xadd = uniform_grid([r1+r, 2*Rx], [-Ry, Ry], h=0.4) + \
       uniform_grid([0, r1], [Ry0, Ry], h=0.4) + \
       uniform_grid([0, r1], [-Ry, -Ry0], h=0.4)
Xdna = uniform_grid([r0-r+.01, r1+r-.01], [-l0/2-r, l0/2+r], h=0.2)
scatter(X, c="b")
scatter(Xadd, c="r")
scatter(Xdna, c="g")

Fexpadd = [Fimpl(x) for x in Xadd]
X += Xadd
Fexp += Fexpadd
Fimp = [Fimpl(x0) for x0 in X]

X += Xdna
Fdna = [[0.,0.] for x in Xdna]
Fexp += Fdna
Fimp += Fdna
#y = 1e-12*np.array(Fexp)/np.array(Fimp) - 1
Fexp = np.array(Fexp)
Fimp = np.array(Fimp)

x = np.array([t[0] for t in X])
y = np.array([t[1] for t in X])

# overwrite explicit colculation far from pore
yfar = abs(y) > Ry0
Fexp[yfar,:] = Fimp[yfar,:]

# duplicate array
notx0 = x>0.
def mirror(z, sign):
    return np.concatenate([z, sign*z[notx0]])

x = mirror(x,-1)
y = mirror(y,1)

trid = dln.Triangulation(x, y)
tri = mtri.Triangulation(x, y)

plt.figure()
plt.triplot(tri, '-')
Fi = [None]*2
mesh = nanopores.RectangleMesh([-2*Rx, -Ry], [2*Rx, Ry], 60, 150)

for i in (0, 1):
    z = Fexp[:,i]
    z = mirror(z, -1 if i==0 else 1)
    
    # interpolate data
    interp = trid.nn_interpolator(z)
    interps = lambda x: interp([x[0]], [x[1]])
    Fi[i] = S1(interps, mesh)
    #dolfin.plot(Fi[i])

V = dolfin.VectorFunctionSpace(mesh, "CG", 1)
FFexp = dolfin.Function(V)
dolfin.assign(FFexp, [Fi[0], Fi[1]])

#dolfin.plot(FFexp)    
poly=Howorka.polygon()
porestreamlines(polygon=poly, rx=4., ry=7., F=FFexp)

#dolfin.interactive()    
plt.show()

#exit()
##zo = interp(xi,yi)
#
## define grid.
#xi = np.linspace(0., Rx, 40)
#yi = np.linspace(-Ry, Ry, 200)
## grid the data.
#zi = griddata(x0, x1, z, xi, yi, interp='nn')
#print list(zi)
## contour the gridded data, plotting dots at the nonuniform data points.
#CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
#CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,
#                  vmax=abs(zi).max(), vmin=-abs(zi).max())
#plt.colorbar()  # draw colorbar
## plot data points.
#plt.scatter(x0, x1, marker='o', c='b', s=5, zorder=10)
##plt.xlim(-2, 2)
##plt.ylim(-2, 2)
#plt.title('griddata test (%d points)' % 114)
#plt.show()

"""

# visualization
alpha = [interpoland(x, [f[i] for f in y], eps=eps) for i in (0,1)]
def FF(x):
    f0 = (1. + alpha[0](x))*Fimpl(x)[0]
    f1 = (1. + alpha[1](x))*Fimpl(x)[1]
    return np.array([f0, 0, f1])
def FF2D(x):
    f0 = (1. + alpha[0](x))*Fimpl(x)[0]
    f1 = (1. + alpha[1](x))*Fimpl(x)[1]
    return np.array([f0, f1])    
    
print FF2D(np.array([0,0]))

#F0 = S1(lambda x: (1. + alpha[0](x))*Fimpl(x)[0], Rx, Ry)
F0 = S1(alpha[0], Rx, Ry)
dolfin.plot(F0)
F2 = S1(alpha[1], Rx, Ry)
dolfin.plot(F2)
#
#FF_S1_2 = S1(FF2D, Rx, Ry, dim=2)
#FF_S1_3 = S1(FF, Rx, Ry, dim=3)
#dolfin.plot(FF_S1_2)

# plot streamlines
poly=Howorka.polygon()
FF_S1_3 = S1(FF, Rx, Ry, dim=3)
porestreamlines(polygon=None, rx=Rx, ry=Ry, F=FF_S1_3)
nanopores.showplots()
dolfin.interactive()
"""
"""
# visualization
Flambdas = [interpoland(x, [f[i] for f in F], eps=eps) for i in (0,1,2)]
def FF2D(x):
    return np.array([Flambdas[0](x), Flambdas[2](x)])
def FF(x):
    return 1e-12*np.array([fi(x) for fi in Flambdas])
print FF2D(np.array([0,0]))

F0 = S1(Flambdas[0], Rx, Ry)
dolfin.plot(F0)
F2 = S1(Flambdas[2], Rx, Ry)
dolfin.plot(F2)

FF_S1_2 = S1(FF2D, Rx, Ry, dim=2)
FF_S1_3 = S1(FF, Rx, Ry, dim=3)
dolfin.plot(FF_S1_2)

# plot streamlines
poly=Howorka.polygon()
porestreamlines(polygon=None, rx=Rx, ry=Ry, F=FF_S1_3)
nanopores.showplots()
dolfin.interactive()
"""