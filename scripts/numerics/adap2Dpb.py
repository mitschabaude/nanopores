from nanopores import *
from nanopores.geometries.curved import Circle
from dolfin import *
from mysolve import adaptive_pb

geo_name = "H_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 1.,
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
cheapest = False,
ref = 6.08430894614e+14, #2.66339790473e+12, 
adaptq = True,
)

geo_params = dict(
x0 = [0.,0.,z0],
rMolecule = 0.5*nm,
moleculeblayer = False,
membraneblayer = False,
)

phys_params = dict(
Membraneqs = -0.0,
Qmol = -1.*qq,
bulkcon = 3e2,
dnaqsdamp = 1.,
bV = bV,
exactMqv = not adaptq,
adaptMqv = adaptq,
)

meshgen_dict = generate_mesh(h, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)

# define circle for molecule
molec = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
# this causes geo to automatically snap boundaries when adapting
geo.curved = dict(moleculeb = molec.snap)

phys = Physics("howorka", geo, **phys_params)
plot(geo.boundaries)
plot(geo.subdomains)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
#StokesProblemAxisymEqualOrder.beta = 1.0 #1e-18
PNPSAxisym.tolnewton = 1e-2

pb = adaptive_pb(geo, phys, cyl=True, frac=.5, Nmax=Nmax, Fpbref=ref, cheapest=cheapest)

print("hmin [nm]: ", geo.mesh.hmin()/nm)
#print phys
plot(geo.boundaries)
#pb.visualize()

pb.estimators["err ref"].plot(rate=-1.)
if not cheapest:
    pb.estimators["rep"].plot(fig=False)
    pb.estimators["err"].plot(fig=False)

pb.estimators["goal"].plot()
if not cheapest:
    pb.estimators["goal ex"].plot(fig=False)
pb.estimators["goal ref"].plot(fig=False)
#saveplots("adap2Dpb")
showplots()
