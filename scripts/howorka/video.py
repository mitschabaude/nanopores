"""plot/save video of molecule convection on PNPS force field"""

import numpy, dolfin, os
from matplotlib import pyplot
import matplotlib.tri as mtri
import matplotlib
import nanopores
from nanopores.physics.convdiff import ConvectionDiffusion

nanopores.add_params(
    log = True,
    video = False,
    levels = 7,
    steps = 100, # timesteps per level for logarithmic time plot
    contourstep = 0.1, # recommended: 0.1 - 1., determines plot speed/smoothness
)

# import force field, mesh etc.
from .forcefield import (function_from_lambda, mesh, geo, phys, Fel, Fdrag,
                        v2d, Ry, Rx)

# video directories
TMPDIR = "/tmp/video/"
VIDDIR = os.path.expanduser("~") + "/presentations/nanopores/"
if video:
    if not os.path.exists(TMPDIR):
        os.makedirs(TMPDIR)
    else:
        for f in os.listdir(TMPDIR):
            os.remove(TMPDIR + f)

# initial condition
N = 10. # number of molecules to diffuse
r = 1. # radius of spherical region where molecules start [nm]
Vol = dolfin.pi*4./3.*r**3 # volume of region [nm**3]
c0 = N/Vol # concentration [1/nm**3]
#c0 = N # concentration [1/Vol]
x0 = numpy.array([0., 8.]) # position of region    
u0f = lambda x: (c0 if sum((x-x0)**2) < r**2 else 0.) # function

u0 = function_from_lambda(u0f)

def convect(geo, phys, Fel, Fdrag, u0, log=False):
    if log:
        frac = 1./steps
        t = 1e-9    
        dt = t*frac
    else:
        t = 1e-7
        dt = 1e-9 
    bc = {} #dict(upperb=dolfin.Constant(0.), lowerb=dolfin.Constant(0.))
    F = dolfin.Constant(1.)*(Fel + dolfin.Constant(3.)*Fdrag)
    pde = ConvectionDiffusion(geo, phys, dt=dt, F=F, u0=u0, bc=bc, cyl=True)
    pde.add_functionals([current])
    yield pde
    if log:
        pde.timerange = nanopores.logtimerange(t, levels=levels, frac=frac,
                                               change_dt=pde.change_dt)
    for t_ in pde.timesteps(t=t):
        pde.record_functionals()
        yield pde
    if log:
        pde.plot_functionals("semilogx")
    else:
        pde.plot_functionals("plot")
    pyplot.ylabel("current through pore [1/ms]")
    
def current(U, geo):
    u, = U
    r2pi = dolfin.Expression("2*pi*x[0]")
    phys = geo.physics
    grad = phys.grad
    D = geo.pwconst("Dtarget")
    kT = dolfin.Constant(phys.kT)
    F = Fel + dolfin.Constant(3.)*Fdrag
    # current density [1/nm**3]*[nm/ns] = [1/(ns*nm**2)]
    j = -D*grad(u) + D/kT*F*u
    #lscale = Constant(phys.lscale)
    L = dolfin.Constant(9.) # pore length
    # current in 1/ns
    J = -j[1]/L *r2pi*geo.dx("pore")
    # current in 1/ms
    J = 1e6*J
    return dict(J=J)

# mesh plot preprocessing
# extract x and y coordinates of nodes
x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
triangles = mesh.cells()
# Create triangulation.
#triang = mtri.Triangulation(x, y, triangles)

# duplicate array
notx0 = x>0.
x2 = -x[notx0]
y2 = y[notx0]
xx = numpy.concatenate([x, x2])
yy = numpy.concatenate([y, y2])
N0 = x.shape[0]
triangles2 = numpy.array(triangles)
# TODO: slow
for i, j in enumerate(numpy.where(notx0)[0]):
    triangles2[triangles2==j] = i+N0
tt = numpy.concatenate([triangles, triangles2])
triang = mtri.Triangulation(xx, yy, tt)
zz = numpy.zeros(xx.shape)

def function2values(u):
    z = u.vector()[v2d]
    zz[:N0] = z
    zz[N0:] = z[notx0]
    return zz

# uneven bounds changes the colormapping:
print(c0)
contours = 10**numpy.arange(-8, numpy.log10(c0), contourstep)
contours = numpy.concatenate((numpy.array([0.]), contours)) #numpy.array([0., 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.00]) # contour ranges
#contours = numpy.array([0., 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1])
#contours = numpy.array([0., 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1])
norm = matplotlib.colors.BoundaryNorm(boundaries=contours, ncolors=256)
norm2 = matplotlib.colors.SymLogNorm(linthresh=1e-7, linscale=1.0, vmin=-0., vmax=c0)

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}\cdot10^{{{}}}$'.format(a,b)
formt = matplotlib.ticker.FuncFormatter(fmt)
#numpy.arange(-0.1,1.1,0.1)

# video!!
pyplot.ion()
#fig = pyplot.figure()
#ax = fig.add_subplot(111)
fig, ax = pyplot.subplots(figsize=(12,10))
#ax.set_aspect('equal')
ax.set_ylim([-Ry, Ry])
ax.set_xlim([-Rx, Rx])
ax.set_title("target molecule concentration")
#cax = fig.add_axes([.7, .1, 0.025, .8])
cax, kw = matplotlib.colorbar.make_axes(ax)

isclear = True
i = 0   
for pde in convect(geo, phys, Fel, Fdrag, u0=u0, log=log):
    t = pde.time[-1] if len(pde.time)>0 else 0.
    if not isclear:
        ax.clear()
        cax.clear()
    z = numpy.abs(function2values(pde.solution)) # dirty trick to hide negative values
    #CS = ax.tripcolor(triang, z, cmap='PuBu_r', norm=norm2)
    CS = ax.tricontourf(triang, z, contours, norm=norm2, cmap='PuBu_r')
    ax.set_title("t = %.1e" % t)
    ax.set_ylim([-Ry, Ry])
    ax.set_xlim([-Rx, Rx])
    #CS = ax.tripcolor(triang, z, norm=norm, cmap='PuBu_r')
    #ax.triplot(triang, '-')
    #fig.colorbar(CS, cax=cax, extend="both", orientation="vertical")
    cb = fig.colorbar(CS, cax=cax, extend="both", orientation="vertical", format=formt)
    #cb.set_label(r"molecules per $\rm{nm}^{3}$")
    #pyplot.colorbar(CS)
    fig.canvas.draw()
    if video:
        fig.savefig(TMPDIR+"%03d.png" % i) #, dpi=256)
    isclear = False
    i += 1    
    #pde.visualize()

#pyplot.tricontourf()
#norm=matplotlib.colors.LogNorm()
#pyplot.tricontourf()
if video:
    #os.system("avconv -i %s%%03d.png %svideo.mp4 -y" % (TMPDIR, VIDDIR))
    os.system("avconv -loglevel quiet -i %s%%03d.png %svideo.mp4 -y" % (TMPDIR, VIDDIR))
    for f in os.listdir(TMPDIR):
        os.remove(TMPDIR + f)
#pyplot.show(block=True)
