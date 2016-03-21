from nanopores import *
from nanopores.geometries.curved import Cylinder, Sphere, Circle
from dolfin import *
from mysolve import adaptive_pbpnps, adaptive_pb

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 1.,
h2D = 1.,
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
Nmax2D = 1e4,
frac = .5,
cheapest = True,
ref = 6.08430894614e+14, #2.66339790473e+12
adaptq = True,
ratio = .01,
)
print PARAMS

geo_params = dict(
x0 = [0.,0.,z0],
rMolecule = 0.5*nm,
lcCenter = 0.5,
lcMolecule = 0.5, #025,
)
geo_params2D = dict(
x0 = [0.,0.,z0],
rMolecule = 0.5*nm,
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

Felref = 1.46603151309
Fdragref = -25.8709397291

meshgen_dict = generate_mesh(h, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)

# 2D geometry as well
generate_mesh(h2D, "H_geo", **geo_params2D)
geo2D = geo_from_name("H_geo", **geo_params)
molec2D = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
geo2D.curved = dict(moleculeb = molec2D.snap)

# define sphere for molecule
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
    outermembraneb = side.snap,
    )

phys2D = Physics("pore_molecule", geo2D, **phys_params)
phys = Physics("pore_molecule", geo, **phys_params)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = .7

#PNPProblem.method["iterative"] = False
PNPProblem.method["kparams"]["relative_tolerance"] = 1e-10
PNPProblem.method["kparams"]["absolute_tolerance"] = 1e-6
PNPProblem.method["kparams"]["nonzero_intial_guess"] = False #True
PNPProblem.method["kparams"]["monitor_convergence"] = False #True
PNPProblem.method["iterative"] = False
StokesProblem.method["iterative"] = False #True
#StokesProblemEqualOrder.beta = 1.
LinearPBProblem.method["ks"] = "bicgstab"
LinearPBProblem.method["kparams"]["relative_tolerance"] = 1e-10
LinearPBProblem.method["kparams"]["absolute_tolerance"] = 1e-6
#LinearPBProblem.method["kparams"]["monitor_convergence"] = True
LinearPBProblem.method["kparams"]["nonzero_intial_guess"] = True
#LinearPBProblem.method["iterative"] = False
#set_log_level(100)

PNPS.tolnewton = 1e-4
PNPSAxisym.tolnewton = 1e-4

# solve 2D
print "\n---- SOLVE 2D PROBLEM ----"
pb2D, pnps2D = adaptive_pbpnps(geo2D, phys2D, cyl=True, frac=frac, Nmax=Nmax2D, 
    Felref=Felref, Fdragref=Fdragref, Fpbref=ref, cheapest=cheapest)   
#pb2D = adaptive_pb(geo2D, phys2D, cyl=True, frac=.5, Nmax=Nmax2D,
#    Fpbref=ref, cheapest=cheapest, ratio=ratio)
mesh2D = geo2D.mesh

# 1D visualization
Rz = geo.params["Rz"]
r0 = geo.params["r0"]
plot1D({"phi (2D)": pb2D.solution}, (-Rz, Rz, 101), "y", dim=2, origin=(r0, 0.))
#plot1D({"phi (2D)": pb2D.solution}, (-Rz, Rz, 101), "y", dim=2, origin=(0., 0.))

# solve 3D
print "\n---- SOLVE 3D PROBLEM ----"
pb, pnps = adaptive_pbpnps(geo, phys, frac=frac, Nmax=Nmax,
    Felref=Felref, Fdragref=Fdragref, Fpbref=ref, cheapest=cheapest)
#pb = adaptive_pb(geo, phys, frac=frac, Nmax=Nmax, Fpbref=ref,
#    mesh2D=mesh2D, cheapest=cheapest, ratio=ratio)

print "hmin [nm]: ", geo.mesh.hmin()/nm

# 2D visualization
phi = pb.functions["primal"]
phi2D = pb2D.functions["primal"]
if mesh2D is not None:
    plot_cross(phi, mesh2D, title="pb primal")
    plot(phi2D, title="pb primal 2D")
    #plot_cross(pb.functions["dual"], mesh2D, title="pb dual")
#plot_sliced(geo)

# 1D visualization
'''
Rz = geo.params["Rz"]
r0 = geo.params["r0"]
plot1D({"phi (2D)": phi2D}, (-Rz, Rz, 101), "y", dim=2, origin=(r0, 0.))
plot1D({"phi (3D)": phi},   (-Rz, Rz, 101), "z", dim=3, origin=(0., r0, 0.), newfig=False)
plot1D({"phi (3D)": phi},   (-Rz, Rz, 101), "z", dim=3, origin=(0., -r0, 0.), newfig=False)
plot1D({"phi (3D)": phi},   (-Rz, Rz, 101), "z", dim=3, origin=(-r0, 0., 0.), newfig=False)
plot1D({"phi (3D)": phi},   (-Rz, Rz, 101), "z", dim=3, origin=(r0, 0., 0.), axlabels=("z [nm]", "potential [V]"),
newfig=False)
'''

# convergence plots
pnps.visualize("pore")
pb.estimators["Fel"].plot()
pb.estimators["Fdrag"].plot(fig=False)
pb.estimators["F"].plot(rate=-1., fig=False)

pb.estimators["err ref"].plot(rate=-2./3.)
if not cheapest:
    pb.estimators["rep"].plot(fig=False)
    pb.estimators["err"].plot(fig=False)
pb2D.estimators["err ref"].plot(rate=-1., fig=False)
"""
pb.estimators["goal"].plot()
if not cheapest:
    pb.estimators["goal ex"].plot(fig=False)
pb.estimators["goal ref"].plot(fig=False)
"""
saveplots("adap3D")
interactive()
showplots()
