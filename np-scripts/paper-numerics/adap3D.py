from nanopores import *
from nanopores.geometries.curved import Cylinder, Sphere
from dolfin import *
from mysolve import adaptive_pbpnps

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 7.,
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
)

geo_params = dict(
x0 = [0.,0.,z0],
rMolecule = 0.5*nm,
#lcCenter = 0.1,
#lcMolecule = 0.05
)

phys_params = dict(
Membraneqs = -0.0,
Qmol = -qq,
bulkcon = 3e2,
dnaqsdamp = .1,
bV = bV,
)

meshgen_dict = generate_mesh(h, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
plot(geo.boundaries)
plot_sliced(geo)

# define circle for molecule
molec = Sphere(R=geo.params["rMolecule"], center=geo.params["x0"])
# define cylinders for inner and outer DNA boundary and side boundary
innerdna = Cylinder(R=geo.params["r0"], L=geo.params["l0"])
outerdna = Cylinder(R=geo.params["r1"], L=geo.params["l0"])
side = Cylinder(R=geo.params["R"], L=2.*geo.params["Rz"])
# this causes geo to automatically snap boundaries when adapting
geo.curved = dict(
    moleculeb = molec.snap,
    innerdnab = innerdna.snap,
    outerdnab = outerdna.snap,
    membranednab = outerdna.snap,
    sideb = side.snap,
    outermembraneb = side.snap,)

phys = Physics("pore_molecule", geo, **phys_params)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
#PNPProblem.method["iterative"] = False
PNPProblem.method["kparams"]["relative_tolerance"] = 1e-10
PNPProblem.method["kparams"]["absolute_tolerance"] = 1e-10
PNPProblem.method["kparams"]["monitor_convergence"] = False #True
StokesProblem.method["iterative"] = False #True
StokesProblemEqualOrder.beta = 1. #True
PNPS.tolnewton = 1e-4

pb, pnps = adaptive_pbpnps(geo, phys, frac=.5, Nmax=Nmax, 
    Felref=1.211487, Fdragref=-7.675373, Fpbref=6.523790e+14)

print "hmin [nm]: ", geo.mesh.hmin()/nm
plot(geo.boundaries)
plot_sliced(geo)
pnps.visualize("pore")
pb.estimators["Fel"].plot()
pb.estimators["Fdrag"].plot(fig=False)
pb.estimators["F"].plot(rate=-1., fig=False)

pb.estimators["rep"].plot()
pb.estimators["err ref"].plot(fig=False)
pb.estimators["err"].plot(rate=-1., fig=False)

#pb.estimators["goal"].plot()
#pb.estimators["goal ex"].plot(fig=False)
#pb.estimators["goal ref"].plot(fig=False)
saveplots("adap2D")
showplots()
