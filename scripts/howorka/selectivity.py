"""plot/save video of molecule convection on PNPS force field"""

import numpy, dolfin
from matplotlib import pyplot
import nanopores
from nanopores.physics.convdiff import ConvectionDiffusion

nanopores.add_params(
    log = True,
    levels = 2,
    t = 1e-9,
    z0 = 40.,
    steps = 100, # timesteps per level for logarithmic time plot
    contourstep = 0.1, # recommended: 0.1 - 1., determines plot speed/smoothness
)

# import force field, mesh etc.
from forcefield import function_from_lambda, mesh, geo, phys, Fel, Fdrag, params

# initial condition
N = 10. # number of molecules to diffuse
r = 5. # radius of spherical region where molecules start [nm]
Vol = dolfin.pi*4./3.*r**3 # volume of region [nm**3]
c0 = N/Vol # concentration [1/nm**3]
#c0 = N # concentration [1/Vol]
x0 = numpy.array([0., z0]) # position of region    
u0f = lambda x: (c0 if sum((x-x0)**2) < r**2 else 0.) # function
u0 = function_from_lambda(u0f)

# total concentration
ctot = dolfin.assemble(u0*dolfin.Expression("2*pi*x[0]")*geo.dx())
print "Total concentration:", ctot, "molecules."

def convect(geo, phys, Fel, Fdrag, u0, t=1e-9, log=False):
    frac = 1./steps
    dt = t/steps
    bc = {} #dict(upperb=dolfin.Constant(0.), lowerb=dolfin.Constant(0.))
    F = dolfin.Constant(1.)*(Fel + dolfin.Constant(3.)*Fdrag)
    pde = ConvectionDiffusion(geo, phys, dt=dt, F=F, u0=u0, bc=bc, cyl=True)
    pde.add_functionals([current, selectivity])
    #yield pde
    pde.timerange = nanopores.logtimerange(t, levels=levels, frac=frac,
                                               change_dt=pde.change_dt)
    for t_ in pde.timesteps(t=t):
        pde.record_functionals()
        pde.visualize()
        
    return pde
    
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
    
def selectivity(U, geo):
    u, = U
    r2pi = dolfin.Expression("2*pi*x[0]")
    # total concentration
    urel = u/dolfin.Constant(ctot)
    c = urel *r2pi*geo.dx() 
    ctop = urel *r2pi*geo.dx("bulkfluidtop")
    cbottom = urel *r2pi*geo.dx("bulkfluidbottom")
    cpore = urel *r2pi*geo.dx("pore")
    return dict(c=c, ctop=ctop, cbottom=cbottom, cpore=cpore)

# compute
pde = convect(geo, phys, Fel, Fdrag, u0=u0, t=t, log=log)

# save results
NAME = "howorka2D_selectivity_Q%.0f" % (params["Qmol"],)
results = dict(
    time = pde.time,
    release = pde.functionals["cbottom"].values,
    current = pde.functionals["J"].values,
)
nanopores.save_stuff(NAME, results, params)

pde.plot_functionals("plot", ["cbottom"])
pyplot.ylabel("% release")
pde.plot_functionals("semilogx", ["J"])
pyplot.ylabel("current through pore [1/ms]")
#pyplot.show(block=True)
