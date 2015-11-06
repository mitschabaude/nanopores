from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
import math

# @Benjamin, Gregor TODO:
# -) check permittivity and surface charge of ahem
# -) what biased voltage to use?

geo_params = dict(
    l3 = 30.,
    l4 = 15.,
    R = 40.,
    x0 = [5., 0., 10.], # |x0| > 2.2
    exit_i = None,
)
phys_params = dict(
    bV = .5,
    ahemqs = 0.02,
    rTarget = 0.5*nm,
    bulkcon = 1000,
)
# TODO: discriminate upper/lower side boundary
exit1 = {"upperb"}
exit2 = {"upperb", "lowerb"}

#StokesProblem.method["lusolver"] = "mumps" # doesn't work
#StokesProblem.method["iterative"] = True

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
meshdict = generate_mesh(10., "aHem", **geo_params)

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

ztop = geo.params["l3"]
zbtm = geo.params["zbtm"] - geo.params["l4"]

def avg(u, meas):
    return assemble(u*meas)/assemble(Constant(1.0)*meas)

def exit_times(tau):
    Tmin = tau(xnear)
    Tmax = tau(xfar)
    Tavg = avg(tau, geo.dS("moleculeb"))
    return (Tmin, Tavg, Tmax)

print
print "--- STATISTICS FOR F=0"
etp_noF = LinearPDE(geo, ExitTimeProblem, phys, F=Constant((0.,0.,0.)), exitb=exitb)
etp_noF.solve(verbose=False)
T_noF = exit_times(etp_noF.solution)
print "\nTime [s] to reach bottom from molecule for F=0: (min, avg, max)"
print T_noF

print "\nTime [s] to reach bottom from pore entrance for F=0:"
print etp_noF.solution([0.,0.,-3.])

t = T_noF[1]
dt = t/100

survival1 = TransientLinearPDE(SurvivalProblem, geo, phys, dt=dt, F=Constant((0.,0.,0.)), exitb=exit1)
survival2 = TransientLinearPDE(SurvivalProblem, geo, phys, dt=dt, F=Constant((0.,0.,0.)), exitb=exit2)

def P(z):
    x = [0., 0., z]    
    return survival1.solution(x) - survival2.solution(x)
    
plotter = TimeDependentPlotter(P, [zbtm, ztop, 200], dt)

for t_ in timerange(t, dt):
    survival1.timestep()
    survival2.timestep()
p1 = survival1.solution


survival2.solve(t=t, visualize=True, verbose=False)
p2 = survival2.solution

print
print "After mean time (%s s) to reach bottom from molecule:" %T_noF[1]
for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
    print "Average survival rate in %s: %.3f percent"%(domain,
        100.*assemble(p*geo.dx(domain))/assemble(Constant(1.0)*geo.dx(domain)))

#print "Physics:"
#for item in phys.__dict__.items():
#    print "%s = %s" %item

print
print "--- CALCULATING F from PNPS"
print
    
pde = PNPS(geo, phys)
pde.solve()
#pde.print_results()

(v, cp, cm, u, p) = pde.solutions(deepcopy=True)
F = phys.Feff(v, u)

for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
    print "Average F in %s:"%domain, assemble(F[2]*geo.dx(domain))/assemble(Constant(1.0)*geo.dx(domain))

#VV = VectorFunctionSpace(geo.mesh, "CG", 1)
#F = project(F, VV)

# solve exit time problem


print
print "--- STATISTICS FOR F=F"

etp = LinearPDE(geo, ExitTimeProblem, phys, F=F, exitb=exitb)
etp.solve(verbose=False)

T = exit_times(etp.solution)
print "\nTime [s] to reach bottom from molecule: (min, avg, max)"
print T

print "\nTime [s] to reach bottom from pore entrance:"
print etp.solution([0.,0.,-3.])



#plot(F, title="F")
#etp.visualize("exittime")

# TIMESTEP
t = T[1]
dt = t/100
survival = TransientLinearPDE(SurvivalProblem, geo, phys, dt=dt, F=F, exitb=exitb)
survival.solve(t=t, visualize=True, verbose=False)

p = survival.solution

print
print "After mean time (%s s) to reach bottom from molecule:" %T[1]
for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
    print "Average survival rate in %s: %.3f percent"%(domain,
        100.*assemble(p*geo.dx(domain))/assemble(Constant(1.0)*geo.dx(domain)))
print

