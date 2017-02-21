# (c) 2016 Gregor Mitscha-Baude
"PNPS solvers and visualization for pugh pore"

import numpy as np
import dolfin
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
import nanopores.physics.simplepnps as simplepnps
import nanopores.tools.solvers as solvers
import nanopores.tools.fields as fields
from nanopores.models.pughpoints import tensorgrid as tensorgrid_

default = nano.Params(
geop = nano.Params(
    dim = 3,
    R = pughpore.params["R"],
    H = pughpore.params["H"],
    l0 = pughpore.params["l0"],
    l1 = pughpore.params["l1"],
    l2 = pughpore.params["l2"],
    l3 = pughpore.params["l3"],
    l4 = pughpore.params["l4"],
    diamPore = None, # will override l0,.. if set
    diamDNA = 2.5, # will override l0,.. if diamPore set
    x0 = pughpore.params["x0"],
    rMolecule = pughpore.params["rMolecule"],
    lcMolecule = pughpore.params["lcMolecule"],
    center_z_at_x0 = False,
),
physp = nano.Params(
    Qmol = 5., # charge of trypsin at pH 8.
    bulkcon = 1000.,
    dnaqsdamp = .5,
    bV = -0.1,
    rDPore = .9,
    bulkbc = True,
    Membraneqs = -0.0,
),
solverp = nano.Params(
    h = 1.5,
    frac = 0.2,
    Nmax = 6e5,
    imax = 30,
    tol = 1e-2,
    cheapest = False,
    stokesiter = False, #True
    diffusivity_data = None,
))
defaultp = default.geop | default.physp

class Setup(solvers.Setup):
    default = default

    def init_geo(self, create_geo=True):
        if self.geop.diamPore is not None:
            diamPore = self.geop.diamPore # inner (effective) pore diameter
            diamDNA = self.geop.diamDNA # dna diameter of outer dna layers
            l0 = diamPore + 6.*diamDNA
            l1 = diamPore + 4.*diamDNA
            l2 = diamPore + 2.*diamDNA
            l3 = diamPore
            l4 = l1
            self.geop.update(l0=l0, l1=l1, l2=l2, l3=l3, l4=l4)
        if not create_geo:
            self.geo = None
            return
        if self.geop.dim == 3:
            geo = pughpore.get_geo(self.solverp.h, **self.geop)
            if geo.params["x0"] is not None:
                molec = nano.curved.Sphere(geo.params["rMolecule"],
                                           geo.params["x0"])
                geo.curved = dict(moleculeb = molec.snap)
        if self.geop.dim == 2:
            geo = pughpore.get_geo_cyl(self.solverp.h, **self.geop)
            if geo.params["x0"] is not None:
                molec = nano.curved.Circle(geo.params["rMolecule"],
                                           geo.params["x0"][::2])
                geo.curved = dict(moleculeb = molec.snap)
        self.geo = geo

    def init_phys(self):
        cyl = self.geop.dim == 2
        self.phys = nano.Physics("pore_mol", self.geo, cyl=cyl, **self.physp)

    def prerefine(self, visualize=False):
        return prerefine(self, visualize=visualize)

class SetupNoGeo(Setup):
    def __init__(self, geop=None, physp=None, solverp=None, **params):
        self.init_params(params, geop=geop, physp=physp, solverp=solverp)
        self.init_geo(create_geo=False)
        self.init_phys()

class Plotter(object):
    def __init__(self, setup=None, dim=3):
        if setup is not None:
            self.geo = setup.geo
            self.dim = setup.phys.dim
        else:
            self.dim = dim
        if self.dim == 3:
            if setup is not None:
                R, H = self.geo.params["R"], self.geo.params["H"]
            else:
                R, H = pughpore.params["R"], pughpore.params["H"]
            self.mesh2D = nano.RectangleMesh([-R,-H/2.], [R, H/2.],
                                                 int(4*R), int(2*H))

    def plot(self, u, title="u", **kwargs):
        if self.dim == 3:
            nano.plot_cross(u, self.mesh2D, title=title, key=title, **kwargs)
        elif self.dim == 2:
            dolfin.plot(u, title=title, key=title, **kwargs)

    def plot_vector(self, u, title="u"):
        if self.dim == 3:
            nano.plot_cross_vector(u, self.mesh2D, title=title, key=title)
        elif self.dim == 2:
            dolfin.plot(u, title=title, key=title)

def solve(setup, visualize=False):
    geo, phys, solverp = setup.geo, setup.phys, setup.solverp
    if visualize:
        plotter = Plotter(setup)

    if geo.mesh.num_cells() < solverp.Nmax:
        pb = prerefine(setup, visualize)
    else:
        pb = None

    it = phys.dim==3
    # solve 1D problem for side BCs
    set_sideBCs(phys, setup.geop, setup.physp)

    # if given, use precomputed ion diffusivity
    set_D_from_data(phys, solverp.diffusivity_data)

    pnps = simplepnps.PNPSFixedPointbV(geo, phys, ipicard=solverp.imax,
               verbose=True, tolnewton=solverp.tol, #taylorhood=True,
               stokesiter=(it and solverp.stokesiter), iterative=it,
               cyl=phys.cyl)

    print "Number of cells:", geo.mesh.num_cells()
    print "DOFs:", pnps.dofs()
    dolfin.tic()
    for i in pnps.fixedpoint(ipnp=6):
        if visualize:
            v, cp, cm, u, p = pnps.solutions()
            plotter.plot(v, "potential")
            #plotter.plot_vector(u, "velocity")
    print "CPU time (solve): %.3g s" %(dolfin.toc(),)
    return pb, pnps

def get_forces(setup, pnps):
    forces = pnps.evaluate(setup.phys.CurrentPNPSDetail)
    forces.update(pnps.evaluate(setup.phys.ForcesPNPS))
    return forces

def prerefine(setup, visualize=False):
    geo, phys, p = setup.geo, setup.phys, setup.solverp
    dolfin.tic()
    if setup.geop.x0 is None:
        goal = phys.CurrentPB
    else:
        goal = lambda v: phys.CurrentPB(v) + phys.Fbare(v, phys.dim-1)
    pb = simplepnps.SimpleLinearPBGO(geo, phys, goal=goal, cheapest=p.cheapest)

    for i in pb.adaptive_loop(p.Nmax, p.frac, verbose=True):
        if visualize:
            if phys.dim==3:
                nano.plot_sliced_mesh(geo, title="adapted mesh", key="b",
                                      elevate=-90. if i==1 else 0.)
                #dolfin.plot(geo.submesh("solid"), key="b",
                #            title="adapted solid mesh")
            if phys.dim==2:
                dolfin.plot(geo.boundaries, key="b", title="adapted mesh",
                            scalarbar=False)
    print "CPU time (PB): %.3g s" %(dolfin.toc(),)
    return pb

def set_D_from_data(phys, data):
    if data is not None:
        func, mesh = fields.get_functions(**data)
        D = func["D"]
        D = dolfin.as_matrix(np.diag([D[i] for i in range(phys.dim)]))
        phys.update(Dp=D, Dm=D)

def solve1D(geop, physp):
    geo = pughpore.get_geo1D(lc=.01, **geop)
    phys = nano.Physics("pore", geo, **physp)
    pnp = nano.solve_pde(simplepnps.SimplePNPProblem, geo, phys)
    return geo, pnp

def visualize1D(geo, pnp):
    v, cp, cm = pnp.solutions()
    h = geo.params["H"]
    nano.plot1D({"potential": v}, (-h/2, h/2, 1001),
                "x", dim=1, axlabels=("z [nm]", "potential [V]"))
    nano.plot1D({"c+": cp, "c-":cm},  (-h/2, h/2, 1001),
                "x", dim=1, axlabels=("z [nm]", "concentrations [mol/m^3]"))

class u1D(dolfin.Expression):
    def __init__(self, u, damping=1., **kw):
        self.u = u
        self.damping = damping
        #dolfin.Expression.__init__(self)

    def damp(self, scalar):
        self.damping *= scalar

    def eval(self, value, x):
        dim = x.shape[0]
        value[0] = self.damping*self.u(x[dim-1])

def set_sideBCs(phys, geop, physp):
    geo, pnp = solve1D(geop, physp)
    v, cp, cm = pnp.solutions()
    phys.v0["sideb"] = u1D(v, degree=1)
    phys.cp0["sideb"] = u1D(cp, degree=1)
    phys.cm0["sideb"] = u1D(cm, degree=1)

def join_dicts(list):
    "[{'F':1.0}, {'F':2.0}, ...] --> {'F':[1.0, 2.0, ...]}"
    return {key:[dic[key] for dic in list] for key in list[0]}

# evaluate finite-size model for a number of x positions
@solvers.cache_forcefield("pugh_force", defaultp)
def F_explicit(X, **params):
    _params = dict(defaultp, **params)
    if "x0" in _params: _params.pop("x0")
    values = []
    for x0 in X:
        setup = Setup(x0=x0, **_params)
        pb, pnps = solve(setup, False)
        values.append(get_forces(setup, pnps))
    return join_dicts(values)

# TODO
## evaluate point-size model for a number of z positions
#def F_implicit3D(Z, **params):
#    geo, phys = setup2D(z0=None, **params)
#    pb, pnps = solve2D(geo, phys, **params)
#    values = [pnps.zforces_implicit(z0) for z0 in Z]
#    F, Fel, Fdrag = tuple(zip(*values))
#    return F, Fel, Fdrag
#
## get discrete force fields from point-size model
#def F_field_implicit3D(**params):
#    params["z0"] = None
#    geo, phys = setup2D(**params)
#    pb, pnps = solve2D(geo, phys, **params)
#    (v, cp, cm, u, p) = pnps.solutions()
#    F, Fel, Fdrag = phys.Forces(v, u)
#    return F, Fel, Fdrag

def tensorgrid(nz=30, nr=4, plot=False, eps=5e-2, eps2=1e-1, buf=7., **params):
    setup = SetupNoGeo(**params)
    return tensorgrid_(nz, nr, plot, eps, eps2, buf, **setup.geop)

def polygon(rmem = 20., **params):
    "polygon of pore + membrane for plotting"
    setup = SetupNoGeo(**params)
    params = nano.Params(pughpore.params) | setup.geop

    r = [0.5*params.l3, 0.5*params.l2, 0.5*params.l1, 0.5*params.l0,
         0.5*params.l4, rmem]
    ztop = params.hpore/2.
    zbot = -ztop
    z = [zbot, ztop - params.h2, ztop - params.h1, ztop, zbot + params.h4,
         zbot + params.hmem]
    # indices: [(0,0), (0,1), (1,1), (1,2), ..., (5,5), (5,0)]
    return [(r[i / 2 % 6], z[(i+1) / 2 % 6]) for i in range(12)]


#........................R.............................
#                                                     .
#                                                     .
#              .........l0..........                  .
#              .                   .                  .
#              ._ _______________ _...............    .
#              |D|               |D|     .   .   .    .
#              |D|......l1.......|D|    h1   .   .    .
#              |D|_ ____l2_____ _|D|......   h2  .    .
#              |DDD|_ _______ _|DDD|..........   .    .
#              |DDDDD|       |DDDDD|             .    .
#              |DDDDD|       |DDDDD|             .    .
#       DNA--->|DDDDD|       |DDDDD|           hpore  .
#              |DDDDD|       |DDDDD|             .    .
#              |DDDDD|..l3...|DDDDD|             .    .
#   MEMBRANE   |DDDDD|       |DDDDD|             .    H
#      |       |DDDDD|       |DDDDD|             .    .
#      |       |DDDDD|       |DDDDD|....h4       .    .
#______V_________|DDD|       |DDD|_____.________ .___ .......
#MMMMMMMMMMMMMMMM|DDD|       |DDD|MMMMM.MMMMMMMMM.MMMM.    hmem
#MMMMMMMMMMMMMMMM|DDD|_______|DDD|MMMMM.MMMMMMMMM.MMMM.......
#                .               .                    .
#                .......l4........                    .
#                                                     .
#                                                     .
#                                                     .
#......................................................

if __name__ == "__main__":
    setup = Setup(h=1., Nmax=2e6, dim=3) #, x0=None)
    _, pnps = solve(setup, True)
    print get_forces(setup, pnps)

    plotter = Plotter(setup)
    v, cp, cm, u, p = pnps.solutions()
    plotter.plot_vector(u, "velocity")
    plotter.plot(cm, "cm")
    plotter.plot(cp, "cp")
    plotter.plot(p, "p")
    dolfin.interactive()
