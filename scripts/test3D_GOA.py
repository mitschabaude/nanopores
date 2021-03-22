from nanopores import *
from dolfin import *
parameters['form_compiler']['representation'] = 'quadrature'

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]
z0 = -4.*nm
#x0 = [4.575, 0., -13.585]

geo_params = dict(
#x0 = None,
x0 = [0.,0.,z0],
#x0 = [4.575, 0., -13.585],
#x0 = [0., 0., -8.372],
rMolecule = 0.55*nm,
#moleculeblayer = True,
)

phys_params = dict(
Membraneqs = -0.,
#bV = 0.011,
Qmol = -2.*qq,
bulkcon = 3e2,
#uppermbias = 1.,
lowermbias = -.02,
dnaqsdamp = 0.5,
)

phys_params_PB = dict(
Membraneqs = -0.,
bV = 0.01,
Qmol = -0.*qq,
bulkcon = 3e2,
#uppermbias = 1.,
#lowermbias = -.02,
dnaqsdamp = 0.,
)

t = Timer("Mesh Generation")
meshgen_dict = generate_mesh(7., geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)


#mesh = geo.submesh("fluid")
#geo = geo_from_name(geo_name, mesh=mesh, **geo_params)
#geo_from_subdomains(mesh, "nanopores.geometries.%s.subdomains" %geo.params["name"], **geo.params)

'''
#plot(geo.boundaries)
#interactive()
for subd in geo._physical_domain:
    #subd = "membrane"
    submesh = geo.submesh(subd)
    geo_sub = geo_from_subdomains(submesh,
                "nanopores.geometries.%s.subdomains" %geo.params["name"], **geo.params)
    #plot(geo_sub.boundaries, title=("boundaries on %s" %subd), elevate=-3e1)
    plot(geo_sub.subdomains, title=("boundaries on %s" %subd), elevate=-3e1)
    #plot(submesh, title=("initial mesh on %s" %subd), wireframe=True, elevate=-3e1)
interactive()
exit()
'''

#plot(geo.pwconst("initial_ions"))
#plot_on_sub(geo.pwconst("permittivity"), geo, "solid")
#interactive()
#exit()

if geo_params["x0"] is None:
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"moleculeb":set()})
    geo.import_synonymes(synonymes)

print("CPU Time (mesh generation):",t.stop())
print("hmin:", geo.mesh.hmin())

IllposedNonlinearSolver.newtondamp = 1.
#IllposedLinearSolver.stab = 1e12
#StokesProblemEqualOrder.beta = 1.0

PNPS.imax = 50
PNPS.tolnewton = 1e-0

#PNPProblem.method["iterative"] = False
#PNPProblem.method["kparams"]["monitor_convergence"] = True
PNPProblem.method["kparams"]["maximum_iterations"] = 1000
#PNPProblem.method["kparams"]["error_on_nonconvergence"] = False
PNPProblem.method["kparams"]["absolute_tolerance"] = PNPS.tolnewton*1e-3
PNPProblem.method["kparams"]["relative_tolerance"] = 1e-8

StokesProblemEqualOrder.method["iterative"] = False
StokesProblemEqualOrder.method["kparams"]["absolute_tolerance"] = PNPS.tolnewton*1e-3
StokesProblemEqualOrder.method["kparams"]["monitor_convergence"] = True
#StokesProblemEqualOrder.method["kp"] = "hypre_amg"

#LinearPBProblem.method["kparams"]["monitor_convergence"] = True
PoissonProblem.method["kparams"]["relative_tolerance"] = 1e-6
LinearPBProblem.method["kparams"]["relative_tolerance"] = 1e-6
LinearPBProblem.method["ks"] = "bicgstab"
#StokesProblem.method["iterative"] = False
#set_log_level(PROGRESS)

#poisson = Poisson(geo, phys=phys)
#poisson.solve()
#v1 = poisson.solution

t = Timer('adaptive PB')
physPB = Physics("pore_molecule", geo, **phys_params_PB)
pb0 = LinearPB(geo, physPB)
pb0.maxcells = 10e4
pb0.marking_fraction = 0.25
pb0.solve(refinement=True)

geo = pb0.geo
phys = Physics("pore_molecule", geo, **phys_params)
print(phys.charge)

goal = (lambda v : phys.Fbare(v, 2) + phys.CurrentPB(v)) if geo.parameter("x0") else (lambda v : phys.CurrentPB(v))
pb = LinearPBGoalOriented(geo, phys, goal=goal)
pb.maxcells = 50e4
pb.marking_fraction = 0.25
pb.solve(refinement=True)

print("CPU Time (PB with adaptivity):",t.stop())
try:
    pb.estimators["err"].plot(rate=-2./3.)
except: pass
#pb.visualize("solid")
geo = pb.geo
v0 = pb.solution

t = Timer('PNPS')

#phys.bV = phys_params["bV"]
pnps = PNPS(geo, phys, v0=v0)
#pnps = PNPS(geo, phys)
pnps.solvers.pop("Stokes")
pnps.solve(visualize=False, print_functionals=True)
print("CPU Time (PNPS):",t.stop())

pnps.print_results()

#for est in pnps.estimators.values():
#    est.plot()

(v,cp,cm,u,p) = pnps.solutions(deepcopy=True)
#plot_on_sub(v, geo, "fluid", title="v")
#plot_on_sub(v0, geo, "fluid", title="v0")
#plot_on_sub(v1, geo, "fluid", title="v1")
#interactive()
l0 = 10*nm
I = pnps.get_functionals()["Javgbtm"]
V = v([0.0, 0.0, -l0]) - v([0.0, 0.0, l0])

print("I (current through pore center):",I,"[pA]")
print("V (transmembrane potential):",V,"[V]")
print("conductance I/V:",I/V,"[pS]")
#print v0([0.0, 0.0, -l0]) - v0([0.0, 0.0, l0])
#print v1([0.0, 0.0, -l0]) - v1([0.0, 0.0, l0])

l1 = 10*nm
cmdiff = cm([0*nm, 0., -l1]) - cm([0*nm, 0., l1])
print("cm difference across membrane:",cmdiff,"[mol/m**3]")

print("u(0) = ",u(0.,0.,0.))

if geo.params["x0"] is None:
    r = geo.params["rMolecule"]
    R = geo.params["r0"]
    if "x0" in globals():
        E0 = (v([x0[0], x0[1], x0[2]-r])-v([x0[0], x0[1], x0[2]+r]))/(2.*r)*phys.lscale
        Fel0 = 1e12*phys_params["Qmol"]*E0
        Fdrag0 = 1e12*6*pi*eta*r*u(x0)[2]*exp(2.*r/(R-r))/phys.lscale
    else:
        Fel0 = 1e12*phys_params["Qmol"]*(v([0.,0.,z0-r])-v([0.,0.,z0+r]))/(2.*r)*phys.lscale
        Fdrag0 = 1e12*6*pi*eta*r*u([0.,0.,z0])[2]*exp(2.*r/(R-r))/phys.lscale
    print("Fbare (theory) [pN]:", Fel0)
    print("Fdrag (theory) [pN]:", Fdrag0)
    print("F [pN]:", Fdrag0 + Fel0)
else:
    fs = pnps.get_functionals()
    Fdrag = fs["Fp2"] + fs["Fshear2"]
    Fel = fs["Fbarevol2"]
    print("phys.Fbare [pN]:",1e12*phys.lscale**(-3)*assemble(phys.Fbare(v, 2)))
    print("Fbare [pN]:", Fel)
    print("Fdrag [pN]:", Fdrag)
    print("F [pN]:", Fdrag + Fel)


#list_timings()
pnps.visualize("solid")
