from importlib import import_module
from nanopores import *
from nanopores.geometries.curved import Cylinder, Sphere, Circle
from .mysolve import pbpnps
import nanopores.tools.solvers as solvers
import dolfin

__all__ = ["setup2D", "solve2D", "setup3D", "solve3D"]

geo_name = "H_geo"
geo_name3D = "H_cyl_geo"
phys_name = "howorka"
params_geo = import_module("nanopores.geometries.%s.params_geo" %geo_name)
nm = params_geo.nm
params_geo3D = import_module("nanopores.geometries.%s.params_geo" %geo_name3D)
nm3D = params_geo3D.nm

add_params(
h = .5,
h3D = 2.,
z0 = 2.*nm,
x0 = None,
bV = -0.1,
Qmol = -1.,
Nmax = 1e4,
Nmax3D = 1e4,
frac = 0.5,
frac3D = 0.2,
cheapest = False,
dnaqsdamp = .25,
rMolecule = 0.5,
Rx = 12.,
Ry = 12.,
bulkcon = 3e2,
stokesLU = False, # determines solver behaviour in 3D
lcpore = 0.5, # relative mesh width of pore in 3D
dim = 2,
)

def geo_params(z0, rMolecule, Rx, Ry):
    x0 = None if z0 is None else [0., 0., z0]
    return dict(
x0 = x0,
rMolecule = rMolecule*nm,
moleculeblayer = False,
membraneblayer = False,
Rx = Rx,
Ry = Ry,
)

def geo_params3D(x0, rMolecule, Rx, Ry, lcpore):
    return dict(
x0 = x0,
rMolecule = rMolecule*nm3D,
lcCenter = lcpore,
lcMolecule = lcpore,
Rx = Rx,
Ry = Ry,
R = Rx,
Rz = Ry,
)

def phys_params(bV, Qmol, dnaqsdamp, rMolecule, bulkcon): return dict(
Membraneqs = -0.0,
Qmol = Qmol*qq,
bulkcon = bulkcon,
dnaqsdamp = dnaqsdamp,
bV = bV,
qTarget = Qmol*qq,
rTarget = rMolecule*1e-9
)

# make compatible with setup paradigm and new style physics
default_geop = Params(geo_params3D(x0, rMolecule, Rx, Ry, lcpore))
default_geop["dim"] = dim
default_physp = Params(phys_params(bV, Qmol, dnaqsdamp, rMolecule, bulkcon))
default_physp["Qmol"] = Qmol
default_physp["bulkbc"] = False
default_solverp = Params(
    h = 2.,
    Nmax = 1e4,
    frac = 0.2,
    cheapest = False,
    stokesLU = False,
)
default = dict(geop=default_geop, physp=default_physp, solverp=default_solverp)

class Setup(solvers.Setup):
    default = default

    def init_geo(self):
        if self.geop.dim == 3:
            generate_mesh(self.solverp.h, geo_name3D, **self.geop)
            geo = geo_from_name(geo_name3D, **self.geop)

            # define cylinders for DNA and side boundaries
            innerdna = Cylinder(R=geo.params["r0"], L=geo.params["l0"])
            outerdna = Cylinder(R=geo.params["r1"], L=geo.params["l0"])
            side = Cylinder(R=geo.params["R"], L=2.*geo.params["Rz"])
            geo.curved = dict(
                innerdnab = innerdna.snap,
                outerdnab = outerdna.snap,
                membranednab = outerdna.snap,
                sideb = side.snap,
                outermembraneb = side.snap,
                )
            # define sphere for molecule
            if self.geop.x0 is not None:
                molec = Sphere(R=geo.params["rMolecule"],
                               center=geo.params["x0"])
                geo.curved.update(moleculeb = molec.snap)

        if self.geop.dim == 2:
            generate_mesh(self.solverp.h, geo_name, **self.geop)
            geo = geo_from_name(geo_name, **self.geop)
            #dolfin.plot(geo.boundaries)
            #dolfin.interactive()
            if self.geop.x0 is not None:
                molec = Circle(geo.params["rMolecule"],
                               geo.params["x0"][::2])
                geo.curved = dict(moleculeb = molec.snap)
        # compatibility with new style physics current
        geo.params["lpore"] = geo.params["l0"]
        self.geo = geo


    def init_phys(self):
        cyl = self.geop.dim == 2
        self.phys = Physics("pore_mol", self.geo, cyl=cyl, **self.physp)
        #self.phys = Physics(phys_name, self.geo, **self.physp)

    def prerefine(self, visualize=False):
        return prerefine(self, visualize=visualize)

from nanopores.models.pughpore import prerefine
#def prerefine_old(setup, visualize=False):
#    geo, phys, p = setup.geo, setup.phys, setup.solverp
#    cyl = phys.dim == 2
#    dolfin.tic()
#
#    LinearPB = LinearPBAxisymGoalOriented if cyl else LinearPBGoalOriented
#    z = phys.dim - 1
#    bV = phys.bV
#    phys.update(bV=0.)
#    goal = (lambda v : phys.Fbare(v, z)) if geo.parameter("x0")\
#            else (lambda v : phys.CurrentPB(v))
#    pb = LinearPB(geo, phys, goal=goal)
#    phys.update(bV=bV)
#    pb.maxcells = p.Nmax
#    pb.marking_fraction = p.frac
#    if p.cheapest:
#        pb.estimate = pb.estimate_cheap
#    refined = True
#    i = 0
#
#    print "Number of cells:", pb.geo.mesh.num_cells()
#    while refined:
#        i += 1
#        print "\nSolving PB."
#        pb.single_solve()
#        if visualize:
#            if phys.dim==3:
#                plot_sliced_mesh(geo, title="adapted mesh", key="b",
#                                      elevate=-90. if i==1 else 0.)
#            if phys.dim==2:
#                dolfin.plot(geo.boundaries, key="b", title="adapted mesh",
#                            scalarbar=False)
#        print "\nError estimation."
#        (ind, err) = pb.estimate()
#        print "\nMesh refinement."
#        refined = pb.refine(ind)
#        if not refined:
#            print "Maximal number of cells reached."
#        else:
#            print "New total number of cells:", pb.geo.mesh.num_cells()
#
#    print "CPU time (PB): %.3g s" %(dolfin.toc(),)
#    return pb

# TODO this should depend on the geo_params
def polygon(rMem = 20.):
    "polygon of pore + membrane for plotting"
    r0 = params_geo.r0 # pore radius
    r1 = params_geo.r1 # barrel outer radius
    l0 = 0.5*params_geo.l0 # half pore length
    l1 = 0.5*params_geo.l1 # half membrane thickness
    return [[r0, l0], [r1, l0], [r1, l1], [rMem, l1],
           [rMem,-l1], [r1, -l1], [r1,-l0], [r0,-l0]]

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
PNPSAxisym.tolnewton = 1e-4
PNPS.tolnewton = 1e-4

def setup2D(mesh=None, **params):
    # pass mesh argument to prevent mesh regeneration
    global z0
    global geo # to prevent garbage collection of mesh
    globals().update(params)
    if z0 is not None:
        z0 = round(z0, 4)
    geop = geo_params(z0, rMolecule, Rx, Ry)
    physp = phys_params(bV, Qmol, dnaqsdamp, rMolecule, bulkcon)

    if mesh is None:
        generate_mesh(h, geo_name, **geop)
        geo = geo_from_name(geo_name, **geop)
    else:
        geo = geo_from_name(geo_name, mesh=mesh, **geop)

    if z0 is not None:
        molec = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
        geo.curved = dict(moleculeb = molec.snap)

    phys = Physics(phys_name, geo, **physp)
    #phys.permittivity = {"default": phys.permittivity["water"]}
    return geo, phys

def solve2D(geo, phys, **params):
    globals().update(params)
    return pbpnps(geo, phys, cyl=True, frac=frac, Nmax=Nmax, cheapest=cheapest)

# evaluate finite-size model for a number of z positions
def F_explicit(z, **params):
    import dolfin
    values = []
    for z0 in z:
        geo, phys = setup2D(z0=z0, **params)
        dolfin.plot(geo.boundaries, key="b", title="boundaries")
        pb, pnps = solve2D(geo, phys, **params)
        dolfin.plot(geo.boundaries, key="b", title="boundaries")
        values.append(pnps.zforces())
    F, Fel, Fdrag = tuple(zip(*values))
    return F, Fel, Fdrag

# evaluate point-size model for a number of z positions
def F_implicit(z, **params):
    geo, phys = setup2D(z0=None, **params)
    pb, pnps = solve2D(geo, phys, **params)
    values = [pnps.zforces_implicit(z0) for z0 in z]
    F, Fel, Fdrag = tuple(zip(*values))
    return F, Fel, Fdrag

# get discrete force fields from point-size model
def F_field_implicit(**params):
    params["z0"] = None
    geo, phys = setup2D(**params)
    pb, pnps = solve2D(geo, phys, **params)
    (v, cp, cm, u, p) = pnps.solutions()
    F, Fel, Fdrag = phys.Forces(v, u)
    return F, Fel, Fdrag

# 3D from here on
def setup3D(mesh=None, **params):
    # pass mesh argument to prevent mesh regeneration
    global x0
    global geo3D # to prevent garbage collection of mesh
    globals().update(params)
    if x0 is not None:
        x0 = [round(t, 4) for t in x0]
    geop = geo_params3D(x0, rMolecule, Rx, Ry, lcpore)
    physp = phys_params(bV, Qmol, dnaqsdamp, rMolecule, bulkcon)

    if mesh is None:
        generate_mesh(h3D, geo_name3D, **geop)
        geo = geo_from_name(geo_name3D, **geop)
    else:
        geo = geo_from_name(geo_name3D, mesh=mesh, **geop)

    # define cylinders for inner and outer DNA boundary and side boundary
    innerdna = Cylinder(R=geo.params["r0"], L=geo.params["l0"])
    outerdna = Cylinder(R=geo.params["r1"], L=geo.params["l0"])
    side = Cylinder(R=geo.params["R"], L=2.*geo.params["Rz"])
    # this causes geo to automatically snap boundaries when adapting
    geo.curved = dict(
        innerdnab = innerdna.snap,
        outerdnab = outerdna.snap,
        membranednab = outerdna.snap,
        sideb = side.snap,
        outermembraneb = side.snap,
        )
    # define sphere for molecule
    if x0 is not None:
        molec = Sphere(R=geo.params["rMolecule"], center=geo.params["x0"])
        geo.curved.update(moleculeb = molec.snap)

    phys = Physics(phys_name, geo, **physp)
    #phys.permittivity = {"default": phys.permittivity["water"]}
    return geo, phys

def solve3D(geo, phys, **params):
    globals().update(params)
    if not stokesLU:
        StokesProblem.method["iterative"] = True
    PNPProblem.method["kparams"]["relative_tolerance"] = 1e-10
    PNPProblem.method["kparams"]["absolute_tolerance"] = 1e-6
    PNPProblem.method["kparams"]["nonzero_intial_guess"] = False
    #StokesProblemEqualOrder.beta = 1.
    StokesProblem.method["kparams"].update(
        monitor_convergence = False,
        relative_tolerance = 1e-10,
        absolute_tolerance = 1e-5,
        maximum_iterations = 2000,
        nonzero_initial_guess = True,
        )
    PNPS.tolnewton = 1e-3
    PNPS.alwaysstokes = True
    return pbpnps(geo, phys, frac=frac3D, Nmax=Nmax3D, cheapest=cheapest)

# evaluate finite-size model for a number of z positions
def F_explicit3D(x, **params):
    import dolfin
    values = []
    for x0 in x:
        geo, phys = setup3D(x0=x0, **params)
        dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
        pb, pnps = solve3D(geo, phys, **params)
        dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
        values.append(pnps.forces())
    F, Fel, Fdrag = tuple(zip(*values))
    return F, Fel, Fdrag

# evaluate point-size model for a number of z positions
def F_implicit3D(z, **params):
    geo, phys = setup2D(z0=None, **params)
    pb, pnps = solve2D(geo, phys, **params)
    values = [pnps.zforces_implicit(z0) for z0 in z]
    F, Fel, Fdrag = tuple(zip(*values))
    return F, Fel, Fdrag

# get discrete force fields from point-size model
def F_field_implicit3D(**params):
    params["z0"] = None
    geo, phys = setup2D(**params)
    pb, pnps = solve2D(geo, phys, **params)
    (v, cp, cm, u, p) = pnps.solutions()
    F, Fel, Fdrag = phys.Forces(v, u)
    return F, Fel, Fdrag