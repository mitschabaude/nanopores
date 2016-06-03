from nanopores import *
from nanopores.geometries.curved import Cylinder, Sphere, Circle
from dolfin import *
from mysolve import * #adaptive_pbpnps, adaptive_pb
from nanopores.physics.simplepnps import *

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 1.,
h2D = .2,
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
Nmax2D = 1e4,
frac = .2,
cheapest = False,
adaptq = True,
ratio = .01,
stokesLU = False,
preadapt = False,
)

geo_params = dict(
x0 = [0.,0.,z0],
rMolecule = 0.5*nm,
lcCenter = 0.05 if preadapt else 0.5,
lcMolecule = 0.05 if preadapt else 0.5,
#moleculeblayer = True,
)
geo_params2D = dict(
x0 = [0.,0.,z0],
rMolecule = 0.5*nm,
)

phys_params = dict(
Membraneqs = -0.0,
Qmol = -1.*qq,
bulkcon = 3e2,
dnaqsdamp = .25,
bV = bV,
exactMqv = not adaptq,
adaptMqv = adaptq,
)

set_log_level(30)
ref = load_Fref()

# 2D geometry
generate_mesh(h2D, "H_geo", **geo_params2D)
geo2D = geo_from_name("H_geo", **geo_params)
molec2D = Circle(R=geo2D.params["rMolecule"], center=geo2D.params["x0"][::2])
geo2D.curved = dict(moleculeb = molec2D.snap)
mesh2D = geo2D.mesh
"""
phys2D = Physics("howorka", geo2D, **phys_params)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
PNPSAxisym.tolnewton = 1e-4

# solve 2D
print "\n---- SOLVE 2D PROBLEM ----"
pb2D, pnps2D = adaptive_pbpnps(geo2D, phys2D, cyl=True, frac=frac, Nmax=Nmax2D, 
    cheapest=cheapest, **ref)
#pb2D = adaptive_pb(geo2D, phys2D, cyl=True, frac=.5, Nmax=Nmax2D,
#    Fpbref=ref, cheapest=cheapest, ratio=ratio)
mesh2D = geo2D.mesh

# define Expression for 3D Stokes initial guess
class w3D(dolfin.Expression):
    def __init__(self, u, p):
        self.u = u
        self.p = p
        dolfin.Expression.__init__(self)
    def eval(self, value, x):
        r = sqrt(x[0]**2 + x[1]**2)
        ux = self.u([r, x[2]])
        px = self.p([r, x[2]])
        if not r==0.:
            value[0] = ux[0]*x[0]/r
            value[1] = ux[0]*x[1]/r
        else:
            value[0] = 0.
            value[1] = 0.
        value[2] = ux[1]
        value[3] = px
    def value_shape(self):
        return (4,)
u, p = pnps2D.solutions("Stokes")
w0 = w3D(u=u, p=p)

# 1D visualization
#Rz = geo2D.params["Rz"]
#r0 = geo2D.params["r0"]
#plot1D({"phi (2D)": pb2D.solution}, (-Rz, Rz, 101), "y", dim=2, origin=(r0, 0.))
#plot1D({"phi (2D)": pb2D.solution}, (-Rz, Rz, 101), "y", dim=2, origin=(0., 0.))
"""
# 3D geometry
meshgen_dict = generate_mesh(h, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)

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

phys = Physics("howorka", geo, **phys_params)

IllposedLinearSolver.stab = 1e9
IllposedNonlinearSolver.newtondamp = 1.

#PNPProblem.method["iterative"] = False
PNPProblem.method["kparams"]["relative_tolerance"] = 1e-10
PNPProblem.method["kparams"]["absolute_tolerance"] = 1e-6
PNPProblem.method["kparams"]["nonzero_intial_guess"] = False #True
PNPProblem.method["kparams"]["monitor_convergence"] = False #True
PNPProblem.method["iterative"] = True #False

if not stokesLU:
    StokesProblem.method["iterative"] = True
#StokesProblemEqualOrder.beta = 1.
StokesProblem.method["kparams"].update(
    monitor_convergence = False,
    relative_tolerance = 1e-10,
    absolute_tolerance = 1e-5,
    maximum_iterations = 2000,
    nonzero_initial_guess = True,
    )

LinearPBProblem.method["ks"] = "bicgstab"
LinearPBProblem.method["kparams"]["relative_tolerance"] = 1e-10
LinearPBProblem.method["kparams"]["absolute_tolerance"] = 1e-6
#LinearPBProblem.method["kparams"]["monitor_convergence"] = True
LinearPBProblem.method["kparams"]["nonzero_initial_guess"] = True
#LinearPBProblem.method["iterative"] = False

PNPS.tolnewton = 1e-4
PNPS.alwaysstokes = True
# test
#w = Function(StokesProblemEqualOrder.space(geo.mesh))
#w.interpolate(w0)
#u, p = w.split()
#plot(u)
#plot(p)
#interactive()
#exit()

# solve 3D
print "\n---- SOLVE 3D PROBLEM ----"
pb, pnps = adaptive_pbpnps(geo, phys, frac=frac, Nmax=Nmax,
    cheapest=cheapest, **ref)
"""
pb = adaptive_pb(geo, phys, frac=frac, Nmax=Nmax, Fpbref=ref["Fpbref"],
    mesh2D=mesh2D, cheapest=cheapest, ratio=ratio)
    
pnps = PNPSHybrid(pb.geo, phys_new, inewton=1, verbose=True, nverbose=True) #, v0=pb.solution)
print "\nSolving PNPS."
dofs = pnps.dofs()
print "  Degrees of freedom: %d" % dofs
pnps.single_solve(damp=.9)
"""
"""
fs = pnps.get_functionals()
Fdrag = fs["Fp%d" %z] + fs["Fshear%d" %z]
Fel = fs["Fbarevol%d" %z]
F = Fdrag + Fel
print "Fbare [pN]:", Fel
print "Fdrag [pN]:", Fdrag
print "F     [pN]:", F
"""
print "hmin [nm]: ", geo.mesh.hmin()/nm

# 2D visualization
#phi = pb.functions["primal"]
#phi2D = pb2D.functions["primal"]
(v, cp, cm, u, p) = pnps.solutions()
if mesh2D is not None:
    plot_cross(v, mesh2D, title="potential")
    #plot_cross(cm, mesh2D, title="cm")
    plot_cross(p, mesh2D, title="p")
    plot_cross_vector(u, mesh2D, title="u")
    #plot(phi2D, title="pb primal 2D")
    #plot_cross(pb.functions["dual"], mesh2D, title="pb dual")
plot_sliced(geo)

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
#pnps.visualize("fluid")
pb.estimators["Fel"].plot()
pb.estimators["Fdrag"].plot(fig=False)
pb.estimators["F"].plot(rate=-0.5, fig=False)

pb.estimators["err ref"].plot(rate=-2./3.)
if not cheapest:
    pb.estimators["rep"].plot(fig=False)
    pb.estimators["err"].plot(fig=False)
#pb2D.estimators["err ref"].plot(rate=-1., fig=False)
"""
pb.estimators["goal"].plot()
if not cheapest:
    pb.estimators["goal ex"].plot(fig=False)
pb.estimators["goal ref"].plot(fig=False)
"""
#saveplots("adap3Diterative", meta=PARAMS)
interactive()
showplots()
