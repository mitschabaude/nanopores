import nanopores, dolfin, numpy, os
from nanopores.tools.functions1D import Geometry1D
from nanopores.physics.convdiff import ConvectionDiffusion, ConvectionDiffusionSteady
from nanopores.models import Howorka
from nanopores import Interval
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import force_profiles
from collections import defaultdict
from functools import partial

nanopores.add_params(
    savefig = False,
    steady = False,
    log = True,
    levels = 1,
    t = 1e-8,
    steps = 100, # timesteps per level for logarithmic time plot
)

# create 1D version of pore geometry
a, b = -1000., 1000.
N = 10000
domain = Interval(a, b)
h = 0.5*Howorka.params_geo.l0 # half pore length
pore = Interval(-h, h)
bulkfluidtop = Interval(h, b)

domain.addsubdomains(   
    #fluid = domain,
    pore = pore,
    #bulkfluid = domain - pore
    bulkfluidtop = bulkfluidtop,
    bulkfluidbottom = domain - pore - bulkfluidtop,
)
domain.addboundaries(
    upperb = domain.boundary("right"),
    lowerb = domain.boundary("left"),
)
domain.synonymes = dict(
    bulkfluid = {"bulkfluidtop", "bulkfluidbottom"},
    fluid = {"bulkfluid", "pore"},
)
geo1D = Geometry1D(N=N, domain=domain)
geo = geo1D.geo

phys = nanopores.Physics("howorka", geo=geo, lscale=1e9, rDtargetPore=0.5)
# initial condition
# concentration in 1/nm**3 
c0 = 1.*(phys.mol*phys.nm**3) # = 1 mM = 1*mol/m**3 
u0 = geo.pwconst("c0", dict(bulkfluidtop = c0, default=0.))
# TEST
#geo1D.plot(u0, "x-")
#plt.show()

# total concentration
r0 = Howorka.params_geo.r0
Across = r0**2 * dolfin.pi
ctot = dolfin.assemble(u0*dolfin.Constant(Across)*geo.dx())
print("Total concentration:", ctot, "molecules.")
print("Total concentration:", c0*(b-h)*dolfin.pi, "molecules.")

#from forcefield import geo, phys, Fel, Fdrag, params
def convect(geo, phys, F, u0, t=1e-9, log=False):
    frac = 1./steps
    dt = t/steps
    bc = {} #dict(upperb=dolfin.Constant(c0), lowerb=dolfin.Constant(0.))
    F = dolfin.as_vector([F])
    pde = ConvectionDiffusion(geo, phys, dt=dt, F=F, u0=u0, bc=bc)
    pde.add_functionals([partial(current, F=F), selectivity])
    #yield pde
    pde.timerange = nanopores.logtimerange(t, levels=levels, frac=frac,
                                               change_dt=pde.change_dt)
    for t_ in pde.timesteps(t=t):
        pde.record_functionals()
        #pde.visualize()
    
    return pde
    
def steadystate(geo, phys, F, u0):
    F = dolfin.as_vector([F])
    bc = dict(upperb=dolfin.Constant(c0), lowerb=dolfin.Constant(c0))
    pde = ConvectionDiffusionSteady(geo, phys, F=F, u0=u0, bc=bc)
    pde.add_functionals([partial(current, F=F), selectivity])
    pde.single_solve()
    return pde
    
def current(U, geo, F):
    u = U[0]
    r2pi = dolfin.Constant(Across)
    phys = geo.physics
    grad = phys.grad
    D = geo.pwconst("Dtarget")
    kT = dolfin.Constant(phys.kT)
    # current density [1/nm**3]*[nm/ns] = [1/(ns*nm**2)]
    j = -D*grad(u) + D/kT*F*u
    #lscale = Constant(phys.lscale)
    L = dolfin.Constant(2.*h) # pore length
    # current in 1/ns
    J = -j[0]/L *r2pi*geo.dx("pore")
    # current in 1/ms
    J = 1e6*J
    return dict(J=J)
    
def selectivity(U, geo):
    u = U[0]
    r2pi = dolfin.Constant(Across)
    # urel = % of total concentration
    urel = u/dolfin.Constant(ctot/100)
    c = urel *r2pi*geo.dx() 
    #c = u *r2pi*geo.dx() 
    ctop = urel *r2pi*geo.dx("bulkfluidtop")
    cbottom = urel *r2pi*geo.dx("bulkfluidbottom")
    cpore = urel *r2pi*geo.dx("pore")
    return dict(c=c, ctop=ctop, cbottom=cbottom, cpore=cpore)
    #return dict(c=c)

def gather_currents(name, rMol, DPore=1.):
    currents = defaultdict(list)
    release = defaultdict(list)
    qmols = []
    #figJ = plt.figure("J")
    
    def avg(J):
        return sum(J)/len(J)
    
    for results in force_profiles.Forces(name):
        qmols.append(results["Q"])
        print("\nQ = ", results["Q"], "\n-------------")
        
        for key in "F", "Fi", "Fi2":
            F = results[key]
            F = force_profiles.function_from_lambda(lambda z: 1e-12*F(z))
            F = geo1D.extend_from(F, -10., 10.)
            
            phys = nanopores.Physics("howorka", geo=geo, rTarget=rMol*1e-9, lscale=1e9, rDtargetPore=DPore)
            
            if steady:
                pde = steadystate(geo, phys, F, u0)
                cbottom = pde.functionals["cbottom"].evaluate()
                current = pde.functionals["J"].evaluate()
                print("%s. current: %.3f 1/ms. release (steady): %.3f %%" % (key, current, cbottom))
            else:
                pde = convect(geo, phys, F, u0=u0, t=t, log=log)
                cbottom = pde.functionals["cbottom"].value()
                #current = avg(pde.functionals["J"].values)
                current = pde.functionals["J"].value()
                print("%s. current: %.3f 1/ms. release after %.1g s: %.3f %%" % (key, current, pde.time[-1], cbottom))
            #pde.plot_functionals("semilogx", ["c", "cbottom", "ctop"])
            currents[key].append(current)
            release[key].append(cbottom)
            #force_profiles.plot_function(F, label="Q="+str(results["Q"]))
            #pde.plot_functionals("semilogx", ["J"], fig=figJ)
            #plt.show()
            #if key=="F":
            #    force_profiles.plot_function(u, label="Q="+str(results["Q"]))
        
        #print "Q %s, J %s, Ji %s, Jib %s" % (
        #qmols[-1], currents["F"][-1], currents["Fi"][-1], currents["Fi2"][-1])
    #plt.show()    
    return qmols, currents, release

#c0 = 1.6605 # [mol/m**3] = 1 molecule per (10nm)**3
#names = {0.25: "r025", 0.5: "r05", 0.2: "r02", 0.4: "r04", 0.75: "r075"}
#items = (0.25, "r025"), (0.5, "r05"), (0.75, "r075")

# diffusivities from Paine, Scherr
items = (0.25, "r025", 0.51), (0.5, "r05", 0.17), (0.75, "r075", 0.027)
    
figures = os.path.expanduser("~") + "/papers/pnps-numerics/figures/"    
plot_current = True
plot_release = False
   
for rMol, name, Dpore in items:
    #plt.figure()
    qmols, currents, release = gather_currents(name, rMol, Dpore)
    #plt.legend()
    
    if plot_current:
        fig, ax = plt.subplots(figsize=(5, 4),)
        ax.plot(qmols, currents["F"], "s-g", label="finite-size")
        ax.plot(qmols, currents["Fi"], "s-b", label="point-size")
        ax.plot(qmols, currents["Fi2"], "s--r", label="hybrid-size")
        #ax.set_yscale("log")
        tick_spacing = 1.
        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        #ax.grid(True, which='both')
        #ax.axhline(y=0, color='#cccccc')
        plt.title("r = %s nm" %rMol)
        plt.xlabel("Molecule charge [q]")
        plt.ylabel("Molecule flux [1/ms]")
        plt.legend(loc = "best")
    
    if plot_release:
        fig, ax = plt.subplots(figsize=(5, 4),)
        ax.plot(qmols, release["F"], "s-g", label="finite-size")
        ax.plot(qmols, release["Fi"], "s-b", label="point-size")
        ax.plot(qmols, release["Fi2"], "s--r", label="hybrid-size")
        #ax.set_yscale("log")
        tick_spacing = 1.
        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
        plt.title("r = %s nm" %rMol)
        plt.xlabel("Molecule charge [q]")
        plt.ylabel("% release")
        plt.legend(loc = "best")

    if savefig:
        fig = plt.gcf()
        #fig.set_size_inches(5, 4)
        plt.savefig(figures + "molcurrent_r%.2f.eps" % rMol, bbox_inches='tight')
    
#plt.show()




