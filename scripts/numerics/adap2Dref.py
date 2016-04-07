from nanopores import *
from nanopores.geometries.curved import Circle
from dolfin import *

geo_name = "H_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
)

geo_params = dict(
#x0 = None,
x0 = [0.,0.,z0],
#x0 = [0., 0., -8.372*nm],
rMolecule = 0.5*nm,
#lcMolecule = nm*0.1,
moleculeblayer = False, #True,
membraneblayer = False,
#Rx = 300*nm,
#Ry = 30*nm,
)

phys_params = dict(
Membraneqs = -0.0,
#bV = .1,
Qmol = -1.*qq,
bulkcon = 3e2,
#lowerpbias = .01,
#lowermbias = -.01,
dnaqsdamp = 1.,
bulkconFluo = 10e-3, # bulk concentration of fluorophore [mol/m**3]
hReservoir = 10*nm, # height of cylindrical upper reservoir [nm]
#applylowerqs = True,
couplebVtoQmol = True,
bV0 = 0.01,
)

meshgen_dict = generate_mesh(1., geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)

print geo.params
# define circle for molecule
molec = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
# this causes geo to automatically snap boundaries when adapting
geo.curved = dict(moleculeb = molec.snap)

phys = Physics("howorka", geo, **phys_params)
print phys.charge

plot(geo.boundaries)
plot(geo.subdomains)
print geo
#plot(geo.pwconst("initial_ions"))
#plot(geo.pwconst("permittivity"))
#interactive()
#exit()

if geo.parameter("x0") is None:
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"moleculeb":set()})
    geo.import_synonymes(synonymes)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
#StokesProblemAxisymEqualOrder.beta = 1.0 #1e-18

PNPSAxisym.imax = 50
PNPSAxisym.tolnewton = 1e-2
PNPSAxisym.alwaysstokes = True
#PNPProblemAxisym.k = 2
PNPProblemAxisym.method["iterative"] = False
PNPProblemAxisym.method["kparams"]["monitor_convergence"] = False

phys.bV = 0.
#pb = LinearPBAxisym(geo, phys)
#goal = (lambda v : phys.Fbare(v, 1)) if geo.parameter("x0") else (lambda v : phys.CurrentPB(v))
goal = (lambda v : phys.Fbare(v, 1) + phys.CurrentPB(v)) if geo.parameter("x0") else (lambda v : phys.CurrentPB(v))
#goal = lambda v : phys.CurrentPB(v)
pb = LinearPBAxisymGoalOriented(geo, phys, goal=goal, ref=-3.71676806993e+17)

#pb.imax = 20 # maximal refinements
pb.maxcells = Nmax
pb.marking_fraction = .5
pb.solve(refinement=True)

geo = pb.geo
plot(geo.boundaries)
v0 = pb.solution

#plot_on_sub(v0, geo, "pore", expr=-grad(v0)[1], title="E")
phys.bV = bV
#pnps = PNPSAxisym(geo, phys)
pnps = PNPSAxisym(geo, phys, v0=v0)
pnps.maxcells = Nmax
pnps.marking_fraction = 1.
#pnps.solvers.pop("Stokes")
refine = False
pnps.solve(refinement=refine, print_functionals=True)

# extrapolation of function: does not increase accuracy at all
# problem probably is that functions are zero on part of domain, so extrapolation doesn't work as intended
"""
print "\nExtrapolating solution..."
tic()
upnp = pnps.functions["PNP"]
ustokes = pnps.functions["Stokes"]
mesh = pnps.geo.mesh
V = FunctionSpace(mesh, "CG", 2)
VV = VectorFunctionSpace(mesh, "CG", 2)
Vpnp = MixedFunctionSpace((V, V, V))
Vstokes = VV*V
Upnp = Function(Vpnp)
Ustokes = Function(Vstokes)
Upnp.extrapolate(upnp)
Ustokes.extrapolate(ustokes)
print "Time for extrapolation:", toc(), "s"
# evaluate functionals at extrapolation
print
Jdic = pnps.functionals
for Jstr in sorted(Jdic):
    J = Jdic[Jstr]
    J.form = replace_function_in_form(J.form, upnp, Upnp)
    J.form = replace_function_in_form(J.form, ustokes, Ustokes)
    print ("%s: " %Jstr) + str(J.evaluate())
"""

(v,cp,cm,u,p) = pnps.solutions(deepcopy=True)
I = pnps.get_functionals()["Javgbtm"]
Ry = geo.params["Ry"]
V = v([0.0, -Ry]) - v([0.0, Ry])

for est in pnps.estimators.values():
    if refine: est.plot(rate=-0.5)
    #else: est.plot()

#print "Convergence rates:\n",pnps.estimators["h1"].rates()

#plot_on_sub(v, geo, "pore", expr=-grad(v), title="E")
#interactive()

print
print "I (current through pore center):",I,"[pA]"
print "V (transmembrane potential):",V,"[V]"
print "conductance I/V:",I/V,"[pS]"

l1 = 10*nm
cmdiff = cm([10*nm, -l1]) - cm([10*nm, l1])
print "cm difference across membrane:",cmdiff,"[mol/m**3]"

print "u(0) = ",u(0.,0.)

if geo.params["x0"] is None:
    r = geo.params["rMolecule"]
    R = geo.params["r0"]
    Fel0 = 1e12*phys_params["Qmol"]*(v([0.,z0-r])-v([0.,z0+r]))/(2.*r)
    Fdrag0 = 1e12*6*pi*eta*r*u([0.,z0])[1]*exp(2.*r/(R-r))
    print "Fbare (theory) [pN]:", Fel0
    print "Fdrag (theory) [pN]:", Fdrag0
    print "F [pN]:", Fdrag0 + Fel0
else:
    fs = pnps.get_functionals()
    Fdrag = fs["Fp1"] + fs["Fshear1"]
    Fel = fs["Fbarevol1"]
    print "phys.Fbare [pN]:",1e12*assemble(phys.Fbare(v, 1))
    print "Fbare [pN]:", Fel
    print "Fdrag [pN]:", Fdrag
    print "F [pN]:", Fdrag + Fel

print "hmin [nm]: ", geo.mesh.hmin()/nm
#plot(pnps.geo.mesh)
#interactive()
pnps.visualize()
pb.visualize()
pb.estimators["err"].plot(rate=-1.)
#pb.estimators["err cheap"].plot(fig=False)
pb.estimators["rep"].plot(fig=False)
pb.estimators["err ref"].plot(fig=False)

pb.estimators["goal"].plot()
pb.estimators["goal ex"].plot(fig=False)
pb.estimators["goal ref"].plot(fig=False)
showplots()
