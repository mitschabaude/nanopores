from nanopores import *
from nanopores.geometries.curved import Cylinder, Sphere
from dolfin import *
from mysolve import adaptive_pbpnps, adaptive_pb

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 7.,
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
frac = .5,
cheapest = True,
ref = 2.66438346218e+12,
)

geo_params = dict(
x0 = [0.,0.,z0],
rMolecule = 0.5*nm,
lcCenter = .5, #0.05,
lcMolecule = .5, #0.025,
)
geo_params2D = dict(
x0 = [0.,0.,z0],
rMolecule = 0.5*nm,
)

phys_params = dict(
Membraneqs = -0.0,
Qmol = -1.*qq,
bulkcon = 3e2,
dnaqsdamp = .0,
bV = bV,
adaptMqv = True,
#exactMqv = True,
)

meshgen_dict = generate_mesh(h, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
plot(geo.boundaries)
plot_sliced(geo)

# get 2D mesh as well
generate_mesh(0.2, "H_geo", **geo_params2D)
mesh2D = get_mesh("H_geo")

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
    #outerdnab = outerdna.snap,
    #membranednab = outerdna.snap,
    #sideb = side.snap,
    #outermembraneb = side.snap,
    )

phys = Physics("pore_molecule", geo, **phys_params)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
#PNPProblem.method["iterative"] = False
PNPProblem.method["kparams"]["relative_tolerance"] = 1e-10
PNPProblem.method["kparams"]["absolute_tolerance"] = 1e-10
PNPProblem.method["kparams"]["monitor_convergence"] = False #True
StokesProblem.method["iterative"] = False #True
StokesProblemEqualOrder.beta = 1. #True
LinearPBProblem.method["ks"] = "bicgstab"
LinearPBProblem.method["kparams"]["relative_tolerance"] = 1e-6
PNPS.tolnewton = 1e-4

#pb, pnps = adaptive_pbpnps(geo, phys, frac=frac, Nmax=Nmax, 
#    Felref=1.211487, Fdragref=-7.675373, Fpbref=6.523790e+14)
pb = adaptive_pb(geo, phys, frac=frac, Nmax=Nmax, Fpbref=ref, mesh2D=mesh2D, cheapest=cheapest)

print "hmin [nm]: ", geo.mesh.hmin()/nm
#print phys
plot_sliced(geo)
interactive()
#pnps.visualize("pore")
#pb.estimators["Fel"].plot()
#pb.estimators["Fdrag"].plot(fig=False)
#pb.estimators["F"].plot(rate=-1., fig=False)

pb.estimators["err ref"].plot(rate=-1.)
if not cheapest:
    pb.estimators["rep"].plot(fig=False)
    pb.estimators["err"].plot(fig=False)

pb.estimators["goal"].plot()
if not cheapest:
    pb.estimators["goal ex"].plot(fig=False)
pb.estimators["goal ref"].plot(fig=False)
#saveplots("adap3D")
showplots()
