from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
import math

# @Benjamin, Gregor TODO:
# -) check permittivity and surface charge of ahem
# -) what biased voltage to use?

geo_params = dict(
    l3 = 60.,
    l4 = 10.,
    R = 60.,
    x0 = [5., 0., 10.], # |x0| > 2.2
    exit_i = 1,
)
phys_params = dict(
    bV = .5,
    ahemqs = 0.01,
    rTarget = 0.5*nm,
    bulkcon = 1000.,
)

badexit = {"upperbulkb"}
goodexit = {"exit"}

#StokesProblem.method["lusolver"] = "mumps" # doesn't work
#StokesProblem.method["iterative"] = True
skip_stokes = True
#IllposedNonlinearSolver.newtondamp = 0.8

print
print "--- INPUT VARIABLES:"
print
print "voltage bias: %.4f mV" %(1000.*phys_params["bV"],)
print "a-Hem surface charge: %.4f C/m^2" %(phys_params["ahemqs"],)
print "upper reservoir dimensions: %d x %d x %d nm" %(geo_params["R"], geo_params["R"], geo_params["l3"])
print "molecule position: %d nm above pore" %(geo_params["x0"][2],)
    
print
print "--- MESHING"
print

t = Timer("meshing")
meshdict = generate_mesh(5., "aHem", **geo_params)

print "Mesh generation time:",t.stop()
#print "Mesh file:",meshdict["fid_xml"]
#print "Mesh metadata:"
#for item in meshdict["meta"].items():
#    print "%s = %s" %item
#print 

t = Timer("reading geometry")
geo = geo_from_xml("aHem")

print "Geo generation time:",t.stop()
#print "Geo params:", geo.params
#print "Geo physical domains:", geo._physical_domain
#print "Geo physical boundaries:", geo._physical_boundary

#plot(geo.boundaries)
#plot(geo.submesh("solid"))
#plot(geo.submesh("exittime"))

phys = Physics("pore_molecule", geo, **phys_params)

x0 = geo.params["x0"]
r0 = math.sqrt(sum(x**2 for x in x0))
rnear = r0 - geo.params["rMolecule"]
rfar = r0 + geo.params["rMolecule"]
xnear = map(lambda x: rnear/r0*x, x0)
xfar = map(lambda x: rfar/r0*x, x0)

def avg(u, meas):
    return assemble(u*meas)/assemble(Constant(1.0)*meas)

def exit_times(tau):
    Tmin = tau(xnear)
    Tmax = tau(xfar)
    Tavg = avg(tau, geo.dS("moleculeb"))
    return (Tmin, Tavg, Tmax)
'''
print
print "--- STATISTICS FOR F=0"
etp_noF = LinearPDE(geo, ExitTimeProblem, phys, F=Constant((0.,0.,0.)), exit=badexit)
etp_noF.solve(verbose=False)
T_noF = exit_times(etp_noF.solution)
print "\nTime [s] to reach top from molecule for F=0: (min, avg, max)"
print T_noF

print "\nTime [s] to reach top from pore entrance for F=0:"
print etp_noF.solution([0.,0.,-3.])

Tbot = avg(etp_noF.solution, geo.dx("fluid_bulk_bottom"))
print "\nTime [s] to reach top from bottom reservoir for F=0:"
print Tbot

t = T_noF[1]*.01
dt = t/100

survival = SuccessfulExit(geo, phys, dt=dt, F=Constant((0.,0.,0.)),
    goodexit=goodexit, badexit=badexit)
survival.timerange = logtimerange(t, levels=4, frac=.05, change_dt=survival.change_dt)
    
survival.solve(visualize=True)
p = survival.solution

print
#print "After t = %s s:" %T_noF[1]
for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
    print "Average successful exit rate in %s: %.2f percent"%(domain,
        100.*avg(p, geo.dx(domain)))
print "Average successful exit rate from molecule: %.2f percent" %(100.*avg(p, geo.dS("moleculeb")), )

#print "Physics:"
#for item in phys.__dict__.items():
#    print "%s = %s" %item
'''
print
print "--- CALCULATING F from PNPS"
print
    
pde = PNPS(geo, phys)
if skip_stokes:
    pde.solvers.pop("Stokes")
pde.solve()
#pde.print_results()

(v, cp, cm, u, p) = pde.solutions(deepcopy=True)
F = phys.Feff(v, u)

for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
    print "Average F in %s:"%domain, assemble(F[2]*geo.dx(domain))/assemble(Constant(1.0)*geo.dx(domain))

#VV = VectorFunctionSpace(geo.mesh, "CG", 1)
#F = project(F, VV)
#plot(F, title="F")

print
print "--- STATISTICS FOR F=F"

etp = LinearPDE(geo, ExitTimeProblem, phys, F=F)
etp.solve(verbose=False)
T = exit_times(etp.solution)
print "\nTime [s] to reach bottom from molecule: (min, avg, max)"
print T
print "\nTime [s] to reach bottom from pore entrance:"
print etp.solution([0.,0.,0.])
print "\nTime [s] to reach bottom from bottom reservoir for F=0:"
print avg(etp.solution, geo.dx("fluid_bulk_bottom"))

steadysurv = LinearPDE(geo, SurvivalProblem, phys, F=F, goodexit=goodexit, badexit=badexit)
steadysurv.solve(verbose=False)
steadysurv.visualize("exittime")
psteady = steadysurv.solution

print
print "Steady solution:"
for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
    print "Average successful exit rate in %s: %.2f percent"%(domain,
        100.*avg(psteady, geo.dx(domain)))
print "Average successful exit rate from molecule: %.2f percent"%(
        100.*avg(psteady, geo.dS("moleculeb")),)
print

# TIMESTEP
t = T[1]*.001
dt = t/100

survival = SuccessfulExit(geo, phys, dt=dt, F=F, goodexit=goodexit, badexit=badexit)
p = survival.solution

avgdS = lambda domain: Functional(p/assemble(Constant(1.0)*geo.dS(domain))*geo.dS(domain))
avgdx = lambda domain: Functional(p/assemble(Constant(1.0)*geo.dx(domain))*geo.dx(domain))

survival.functionals = {
    "P(t) from upper reservoir": avgdx("fluid_bulk_top"),
    "P(t) from molecule": avgdS("moleculeb"),
    "P(t) from pore top": avgdx("poretop"),
    "P(t) from pore center": avgdx("porecenter"),
    "P(t) from pore bottom": avgdx("porebottom"),
}
survival.timerange = logtimerange(t, levels=3, frac=.05, change_dt=survival.change_dt)
survival.solve(visualize=True)
survival.plot_functionals("loglog")

print
print "Transient solution after t=%s:" %(survival.time[-1],)
for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
    print "Average successful exit rate in %s: %.2f percent"%(domain,
        100.*avg(p, geo.dx(domain)))
print "Average successful exit rate from molecule: %.2f percent"%(
        100.*avg(p, geo.dS("moleculeb")),)
print

#interactive()

