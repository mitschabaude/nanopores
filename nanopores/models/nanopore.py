# (c) 2017 Gregor Mitscha-Baude
"PNPS solvers and visualization for arbitrary pore"

import numpy as np
import dolfin
import nanopores as nano
import nanopores.physics.simplepnps as simplepnps
import nanopores.tools.solvers as solvers
import nanopores.tools.fields as fields
from nanopores.tools.utilities import tic, toc
from nanopores.geometries.allpores import get_geo

default = nano.Params(
geop = nano.Params(
    geoname = "wei",
    dim = 2,
    x0 = None,
    rMolecule = 0.5,
    receptor = None,
    lcMolecule = 0.1,
    lcCenter = 0.4,
),
physp = nano.Params(
    Qmol = -1.,
    bulkcon = 1000.,
    dnaqsdamp = .5,
    bV = -0.1,
    rDPore = .9,
    Membraneqs = -0.0,
    ahemqs = None,
    ahemuniformqs = False,
    posDTarget = True,
),
solverp = nano.Params(
    h = 1.,
    frac = 0.2,
    Nmax = 2e4,
    imax = 30,
    tol = 1e-2,
    cheapest = False,
    stokesiter = False, #True
    diffusivity_data = None,
    reconstruct = False,
    fluid = "fluid",
    hybrid = True,
))
defaultp = default.geop | default.physp

class Setup(solvers.Setup):
    default = default

    def __init__(self, create_geo=True,
                 geop=None, physp=None, solverp=None, **params):
        self.init_params(params, geop=geop, physp=physp, solverp=solverp)
        self.init_geo(create_geo)
        self.init_phys()

    def init_geo(self, create_geo=True):
        h = self.solverp.h
        re = self.solverp.reconstruct
        self.geo = get_geo(h=h, reconstruct=re,
                           **self.geop) if create_geo else None

    def init_phys(self):
        cyl = self.geop.dim == 2
        if self.geo is None:
            self.physp["rTarget"] = self.geop["rMolecule"]*1e-9
        self.phys = nano.Physics("pore_mol", self.geo, cyl=cyl, **self.physp)

    def prerefine(self, visualize=False):
        return prerefine(self, visualize=visualize)
    
def get_active_params(params):
    setup = Setup(create_geo=False, **params)
    return setup.active_params

class Plotter(object):
    def __init__(self, setup=None, dim=3):
        if setup is not None:
            self.geo = setup.geo
            self.dim = setup.phys.dim
        else:
            self.dim = dim
        if self.dim == 3:
            if setup is not None:
                R = self.geo.params["R"]
                if "Hbot" in self.geo.params:
                    Htop = self.geo.params["Htop"]
                    Hbot = self.geo.params["Hbot"]
                else:
                    H = self.geo.params["H"]
                    Htop = Hbot = H/2.
            else:
                R, Htop, Hbot = 10., 15., 15.
            self.mesh2D = nano.RectangleMesh([-R,-Hbot], [R, Htop],
                                             int(4*R), int(2.*(Htop + Hbot)))

    def plot(self, u, title="u", **kwargs):
        if self.dim == 3:
            nano.plot_cross(u, self.mesh2D, title=title, key=title, **kwargs)
        elif self.dim == 2:
            dolfin.plot(u, title=title, key=title, **kwargs)

    def plot_vector(self, u, title="u", **kwargs):
        if self.dim == 3:
            nano.plot_cross_vector(u, self.mesh2D, title=title, key=title, **kwargs)
        elif self.dim == 2:
            dolfin.plot(u, title=title, key=title, **kwargs)

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

    if solverp.hybrid:
        pnps = simplepnps.PNPSHybrid(geo, phys, ipicard=solverp.imax,
                   verbose=True, nverbose=True, tolnewton=solverp.tol, #taylorhood=True,
                   stokesiter=(it and solverp.stokesiter), iterative=it,
                   cyl=phys.cyl, fluid=solverp.fluid)
    else:
        pnps = simplepnps.PNPSFixedPointbV(geo, phys, ipicard=solverp.imax,
                   verbose=True, tolnewton=solverp.tol, #taylorhood=True,
                   stokesiter=(it and solverp.stokesiter), iterative=it,
                   cyl=phys.cyl, fluid=solverp.fluid)
    #v, cp, cm, u, p = pnps.solutions()
    #plotter.plot(cm, "cm", interactive=True)

    print("Number of cells:", geo.mesh.num_cells())
    print("DOFs:", pnps.dofs())
    tic()
    for i in pnps.fixedpoint(): #ipnp=6):
        if visualize:
            v, cp, cm, u, p = pnps.solutions()
            plotter.plot(v, "potential")
            #plotter.plot(cm, "cm")
            #R, H = tuple(setup.geo.params[s] for s in ["R", "H"])
            #nano.plot1D(dict(cm=cm, cp=cp), dim=2, axis="y",
            #            rng=(-H/2., H/2., 100), origin=[R, 0])
            #dolfin.interactive()
            #nano.showplots()
            #plotter.plot_vector(u, "velocity")
    print("CPU time (solve): %.3g s" %(toc(),))
    return pb, pnps

def get_forces(setup, pnps):
    forces = pnps.evaluate(setup.phys.CurrentPNPSDetail)
    forces.update(pnps.evaluate(setup.phys.ForcesPNPS))
    return forces

def prerefine(setup, visualize=False, debug=False):
    geo, phys, p = setup.geo, setup.phys, setup.solverp
    tic()
    if setup.geop.x0 is None:
        goal = phys.CurrentPB
    else:
        goal = lambda v: phys.CurrentPB(v) + phys.Fbare(v, phys.dim-1)
    if debug:
        plotter = Plotter(setup)
    pb = simplepnps.SimpleLinearPBGO(geo, phys, goal=goal, cheapest=p.cheapest)

    for i in pb.adaptive_loop(p.Nmax, p.frac, verbose=True):
        if visualize:
            if phys.dim==3:
                nano.plot_sliced_mesh(geo, title="adapted mesh", key="b",
                                      elevate=-90. if i==1 else 0.)
            if phys.dim==2:
                dolfin.plot(geo.boundaries, key="b", title="adapted mesh",
                            scalarbar=False)
        if debug:
            u = pb.solution
            plotter.plot(u, "pb")
            #dolfin.interactive()
    print("CPU time (PB): %.3g s" %(toc(),))
    return pb

def set_D_from_data(phys, data):
    if data is not None:
        func, mesh = fields.get_functions(**data)
        D = func["D"]
        #dolfin.plot(D[0], title="Dx", interactive=True)
        D = dolfin.as_matrix(np.diag([D[i] for i in range(phys.dim)]))
        phys.update(Dp=D, Dm=D)

def solve1D(geop, physp):
    geop = dict(geop, dim=1)
    geo = get_geo(h=.01, **geop)
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

# point-size/implicit fields
# TODO: better, GENERAL interface (like cache_forcefield)
#       that can at least control name
def force_pointsize(**params):
    name = "force_pointsize"
    if not fields.exists(name, **params):
        setup = Setup(**params)
        #print setup.geo
        #setup.geo.plot_boundaries(interactive=True)
        _, pnps = solve(setup, True)
        v, cp, cm, u, p = pnps.solutions()
    
        F, Fel, Fdrag = setup.phys.ForceField(v, u, "fluid")
        fields.save_functions(name, setup.active_params,
                              F=F, Fel=Fel, Fdrag=Fdrag)
        fields.update()
    F, = fields.get_functions(name, "F", **params)
    return F

def diffusivity_simple(**params):
    from nanopores.models.diffusion_interpolation import diffusivity_field
    name = "diffusivity_simple"
    if not fields.exists(name, **params):
        setup = Setup(**params)
        dic = diffusivity_field(setup, r=params["rMolecule"],
                                boundary="poresolidb")
        fields.save_functions(name, setup.active_params, **dic)
        fields.update()
    D, = fields.get_functions(name, "D", **params)
    return D

def force_diff(**params):
    # for random walk (with pointsize force field, no current)
    setup = Setup(create_geo=False, **params)
    # DEBUG
    #print "active params", setup.active_params
    #print "inactive params", setup.inactive_params
    F = force_pointsize(**params)
    if setup.phys.posDTarget:
        D = diffusivity_simple(**params)
        name = "diffusivity_div_simple"
        if not fields.exists(name, **params):
            V = D.function_space()
            divD = dolfin.project(dolfin.as_vector([
                      dolfin.grad(D[0])[0], dolfin.grad(D[1])[1]]), V)
            fields.save_functions(name, setup.active_params, divD=divD)
            fields.update()
        divD, = fields.get_functions(name, "divD", **params)
    else:
        D0 = setup.phys.DTargetBulk
        D0a = np.array([D0, D0, D0])
        divDa = np.array([0., 0., 0.])
        D = lambda x: D0a
        divD = lambda x: divDa
    return F, D, divD

# evaluate finite-size model for a number of x positions
@solvers.cache_forcefield("force", defaultp)
def F_explicit(X, **params):
    _params = dict(defaultp, **params)
    if "x0" in _params: _params.pop("x0")
    values = []
    for x0 in X:
        setup = Setup(x0=x0, **_params)
        pb, pnps = solve(setup, False)
        values.append(get_forces(setup, pnps))
    return join_dicts(values)

@solvers.cache_forcefield("IV", defaultp)
def IV(V, **params):
    params["x0"] = None
    for bV, result in nano.collect_dict(V):
        params["bV"] = bV
        setup = Setup(**params)
        pb, pnps = solve(setup, visualize=True)
        result.new = pnps.evaluate(setup.phys.CurrentPNPSDetail)
    return result

# I for different surface charges
@solvers.cache_forcefield("Irho", defaultp)
def Irho(Rho, **params):
    params["x0"] = None
    for rho, result in nano.collect_dict(Rho):
        params["dnaqsdamp"] = rho
        setup = Setup(**params)
        pb, pnps = solve(setup, visualize=True)
        result.new = pnps.evaluate(setup.phys.CurrentPNPSDetail)
    return result

# I for molecule at different z-positions = 1D current profile
@solvers.cache_forcefield("Iz", defaultp)
def Iz(Z, **params):
    for z, result in nano.collect_dict(Z):
        params["x0"] = [0., 0., z]
        setup = Setup(**params)
        #setup.geo.plot_subdomains()
        #setup.geo.plot_boundaries()
        #print setup.geo
        pb, pnps = solve(setup, visualize=True)
        result.new = pnps.evaluate(setup.phys.CurrentPNPSDetail)
    return result

if __name__ == "__main__":
    params = nano.any_params(
        h = 5.,
        Nmax = 4e4,
        dim = 2,
        x0 = [0, 0, -45.],
        rMolecule = 4.,
        lcMolecule = 0.05,
    )
    #ddata = dict(name="Dalphahem", dim=2, Nmax=.4e5, h=1., ahemqsuniform=True, rMolecule=0.11)
    print(params)
    geop = dict(params)
    geop.pop("h")
    setup = Setup(geop=geop, **params)

    #geo, pnp = solve1D(setup.geop, setup.physp)
    #visualize1D(geo, pnp)
    #nano.showplots()

    print(setup.geo.params)
    print(setup.physp)
    print(setup.solverp)
    setup.geo.plot_subdomains()
    #exit()
    #print setup.phys.permittivity
    #dolfin.plot(setup.geo.pwconst("permittivity"), interactive=True)
    _, pnps = solve(setup, True)
    print(get_forces(setup, pnps))

    plotter = Plotter(setup)
    v, cp, cm, u, p = pnps.solutions()
    plotter.plot_vector(u, "velocity")
    plotter.plot(cm, "cm")
    plotter.plot(cp, "cp")
    #plotter.plot(p, "p")
    dolfin.interactive()
