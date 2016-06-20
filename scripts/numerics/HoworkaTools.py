"""do stuff on top of PNPS solve,
like convecting analyte concentration with force field."""

import numpy, dolfin
from matplotlib import pyplot
import matplotlib.tri as mtri
import matplotlib
import nanopores
import Howorka
from nanopores.physics.convectiondiffusion import ConvectionDiffusion

nanopores.add_params(
    log = True,
    save = False,
    **Howorka.PARAMS)
PARAMS.pop("z0")

# save and load implicit force field
FNAME = "howorka2D_implicit"

def save_force_field(**params):
    F, Fel, Fdrag = Howorka.F_field_implicit(**params)
    mesh = Howorka.geo.mesh
    nanopores.save_functions(FNAME, mesh, F=F, Fel=Fel, Fdrag=Fdrag)
    
if save:
    save_force_field(**PARAMS)
    exit()
    
def load_force_field():
    forces, mesh = nanopores.load_vector_functions(FNAME)
    #mesh = forces["F"].function_space().mesh()
    return forces["F"], forces["Fel"], forces["Fdrag"], mesh

F, Fel, Fdrag, mesh = load_force_field()
geo, phys = Howorka.setup2D(mesh=mesh, z0=None, **PARAMS)
submesh = geo.submesh("fluid")

mesh = submesh
geo, phys = Howorka.setup2D(mesh=mesh, z0=None, **PARAMS)
V = dolfin.FunctionSpace(mesh, "CG", 1)
VV = dolfin.VectorFunctionSpace(mesh, "CG", 1)
F = dolfin.interpolate(F, VV)
Fel = dolfin.interpolate(Fel, VV)
Fdrag = dolfin.interpolate(Fdrag, VV)

v2d = dolfin.vertex_to_dof_map(V)
coord = mesh.coordinates() # numpy array of 2-1 arrays

def function_from_values(values):
    u = dolfin.Function(V)
    u.vector()[v2d] = numpy.array(values)
    return u    
def function_from_lambda(f):
    values = [f(x) for x in coord]
    return function_from_values(values)

N = 10. # number of molecules to diffuse
r = 1. # radius of spherical region where molecules start [nm]
Vol = dolfin.pi*4./3.*r**3 # volume of region [nm**3]
c0 = N/Vol # concentration [1/nm**3]
#c0 = N # concentration [1/Vol]
x0 = numpy.array([0., 10.]) # position of region    
u0f = lambda x: (c0 if sum((x-x0)**2) < r**2 else 0.) # function

u0 = function_from_lambda(u0f)

def convect(geo, phys, Fel, Fdrag, u0, log=False):
    if log:
        frac = .05
        t = 1e-9    
        dt = t*frac
        levels = 8
    else:
        t = 1e-7
        dt = 1e-10 
    bc = {} #dict(upperb=dolfin.Constant(0.), lowerb=dolfin.Constant(0.))
    F = dolfin.Constant(10.)*(Fel + dolfin.Constant(3.)*Fdrag)
    pde = ConvectionDiffusion(geo, phys, dt=dt, F=F, u0=u0, bc=bc, cyl=True)
    pde.add_functionals([current])
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
triang = mtri.Triangulation(x, y, triangles)

# uneven bounds changes the colormapping:
print c0
#contours = numpy.exp(numpy.arange(-10, -2))
#contours = numpy.concatenate((numpy.array([0.]), contours)) #numpy.array([0., 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.00]) # contour ranges
contours = numpy.array([0., 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1])
norm = matplotlib.colors.BoundaryNorm(boundaries=contours, ncolors=256)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'$10^{{{}}}$'.format(b)
formt = matplotlib.ticker.FuncFormatter(fmt)
#numpy.arange(-0.1,1.1,0.1)

# video!!
pyplot.ion()
#fig = pyplot.figure()
#ax = fig.add_subplot(111)
fig, ax = pyplot.subplots(figsize=(10,10))
ax.set_aspect('equal')
ax.set_title("target molecule concentration")
#cax = fig.add_axes([.7, .1, 0.025, .8])
cax, kw = matplotlib.colorbar.make_axes(ax)

isclear = True
   
for pde in convect(geo, phys, Fel, Fdrag, u0=u0, log=log):
    t = pde.time[-1]
    if not isclear:
        ax.clear()
        cax.clear()
    z = numpy.abs(pde.solution.vector()[v2d])
    CS = ax.tripcolor(triang, z, norm=norm, cmap='PuBu_r')
    #CS = ax.tricontourf(triang, z, contours, norm=norm, cmap='PuBu_r')
    #ax.triplot(triang, 'k-')
    fig.colorbar(CS, cax=cax, extend="both", orientation="vertical", format=formt)
    #pyplot.colorbar(CS)
    fig.canvas.draw()
    isclear = False    
    #pde.visualize()

#pyplot.tricontourf()
#norm=matplotlib.colors.LogNorm()
#pyplot.tricontourf()
pyplot.show(block=True)
