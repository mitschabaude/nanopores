"""run diffusion equation to determine selectivity of fluophore,
i.e current and release time series for specific molecule."""
from nanopores.tools.fields import cache
import nanopores
import dolfin
import matplotlib.pyplot as plt
from nanopores.physics.convdiff import ConvectionDiffusion
import forcefields

p = nanopores.user_params(
    overwrite = False,
    levels = 1,
    t = 1e-0,
    steps = 100,
    Qmol = -1,
    rMolecule = 0.5,
    implicit = False,
    R = 12.
)

# force field parameters
f_params = dict(
    Qmol = p.Qmol,
    rMolecule = p.rMolecule,
    implicit = p.implicit,
    Ry = p.R,
    Rx = p.R,
    Nmax = 1e5,
)

# parameters for selectivity calculation
sel_params = dict(
    fluocon = 100., # initial concentration [mM] in upper reservoir
    # parameters regarding timestepping
    levels = p.levels, # levels > 1 --> logarithmic time
    t = p.t, # total time of first level
    steps = p.steps, # timesteps per level
)
default = dict(sel_params, **f_params)

def calculate_selectivity(F, geo, phys, fluocon=1, t=1e0, steps=100, levels=1):
    "core functionality of the module"
    
    # concentration in 1/nm**3 (1 M = 0.6 /nm**3)
    c0 = fluocon*(phys.mol*phys.nm**3) 
    u0 = geo.pwconst("c0", dict(bulkfluidtop = c0, default=0.))
    
    # total concentration
    ctot = dolfin.assemble(u0*dolfin.Expression("2*pi*x[0]")*geo.dx())
    phys.ctot = ctot
    print "Total concentration:", ctot, "molecules."
    
    # convect
    phys.F = F
    frac = 1./steps
    dt = t/steps
    bc = {} #dict(upperb=dolfin.Constant(0.), lowerb=dolfin.Constant(0.))
    pde = ConvectionDiffusion(geo, phys, dt=dt, F=F, u0=u0, bc=bc, cyl=True)
    pde.add_functionals([current, concentration])
    pde.timerange = nanopores.logtimerange(t,
       levels=levels, frac=frac, change_dt=pde.change_dt)
    for t_ in pde.timesteps(t=t):
        pde.record_functionals()
        pde.visualize()
    # obtain current, release
    return dict(
        time = pde.time,
        release = pde.functionals["cbottom"].values,
        current = pde.functionals["J"].values)

# functionals    
def current(U, geo):
    u, = U
    r2pi = dolfin.Expression("2*pi*x[0]")
    phys = geo.physics
    grad = phys.grad
    D = geo.pwconst("Dtarget")
    kT = dolfin.Constant(phys.kT)
    F = phys.F
    # current density [1/nm**3]*[nm/ns] = [1/(ns*nm**2)]
    j = -D*grad(u) + D/kT*F*u
    #lscale = Constant(phys.lscale)
    L = dolfin.Constant(9.) # pore length
    # current in 1/ns
    J = -j[1]/L *r2pi*geo.dx("pore")
    # current in 1/ms
    J = 1e6*J
    return dict(J=J)
    
def concentration(U, geo):
    u, = U
    ctot = geo.physics.ctot
    r2pi = dolfin.Expression("2*pi*x[0]")
    # urel = % of total concentration
    urel = u/dolfin.Constant(ctot/100.)
    c = urel *r2pi*geo.dx() 
    ctop = urel *r2pi*geo.dx("bulkfluidtop")
    cbottom = urel *r2pi*geo.dx("bulkfluidbottom")
    cpore = urel *r2pi*geo.dx("pore")
    return dict(c=c, ctop=ctop, cbottom=cbottom, cpore=cpore)

def _diff(dic, keys):
    dic = dic.copy()
    return {k : dic.pop(k) for k in keys}, dic
    
# user interface
@cache("selectivity", default, overwrite=p.overwrite)
def selectivity(params):
    # filter out selectivity params
    sparams, fparams = _diff(params, sel_params.keys())
    
    F, geo, phys = forcefields.F_geo_phys(**fparams)
    result = calculate_selectivity(F, geo, phys, **sparams)
    result["params"] = params
    return result

if __name__ == "__main__":
    import numpy
    results = nanopores.Params(selectivity(**default))
    t = results.time
    J = results.current
    rel = results.release
    params = results.params
    
    plt.figure(0)
    plt.semilogx(t, rel, "x-")
    plt.xlabel("time [s]")
    plt.ylabel("% release")
    plt.title("reservoir size: %.0f nm" % (params["Ry"],))
    plt.ylim(ymin=0.)
    
    def avg(J):
        n = len(J)
        J0 = list(numpy.array(J)[n*0.2:n*0.5])
        return sum(J0)/len(J0)

    plt.figure(1)
    plt.semilogx(t, J, "x-")
    plt.xlabel("time [s]")
    plt.ylabel("current through pore [1/ms]")
    J0 = avg(J)
    plt.plot(t, [J0]*len(t), "k--")
    plt.title("quasi-equilibrium current: %.1f" % J0)
    plt.ylim(ymin=0.)
    
    plt.show()
