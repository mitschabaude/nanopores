from nanopores import *
from nanopores.geometries.curved import Circle
from dolfin import *
from mysolve import * 

geo_name = "H_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 1.,
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
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
dnaqsdamp = .25,
bV = bV,
)

ref = dict(
    Fpbref= 6.08430894614e+14,
    Felref = 1.1424041805,
    Fpref = -3.51633193736,
    Fsref = -4.13093729278,
)
ref = load_Fref()


meshgen_dict = generate_mesh(h, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)

# define circle for molecule
molec = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
# this causes geo to automatically snap boundaries when adapting
geo.curved = dict(moleculeb = molec.snap)

phys = Physics("pore_molecule", geo, **phys_params)
plot(geo.boundaries)
plot(geo.subdomains)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
#StokesProblemAxisymEqualOrder.beta = 1.0 #1e-18
PNPSAxisym.tolnewton = 1e-4

pb, pnps = adaptive_pbpnps(geo, phys, cyl=True, frac=.5, Nmax=Nmax, **ref)

print "hmin [nm]: ", geo.mesh.hmin()/nm
plot(geo.boundaries)
#pnps.visualize()
interactive()
pb.estimators["Fel"].plot()
pb.estimators["Fs"].plot(fig=False)
pb.estimators["Fp"].plot(fig=False)
pb.estimators["F"].plot(rate=-1., fig=False)

pb.estimators["rep"].plot()
pb.estimators["err ref"].plot(fig=False)
pb.estimators["err"].plot(rate=-1., fig=False)

#pb.estimators["goal"].plot()
#pb.estimators["goal ex"].plot(fig=False)
#pb.estimators["goal ref"].plot(fig=False)
#save_Fref(pb, pnps)
#saveplots("adap2D")
showplots()
