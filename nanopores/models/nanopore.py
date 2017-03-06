# (c) 2017 Gregor Mitscha-Baude
"PNPS solvers and visualization for arbitrary pore"

import numpy as np
import dolfin
import nanopores as nano
import nanopores.physics.simplepnps as simplepnps
import nanopores.tools.solvers as solvers
import nanopores.tools.fields as fields
# default geometry
from nanopores.geometries.alphahem import get_geo as get_geo_default

default = nano.Params(
geop = nano.Params(
    dim = 2,
    R = 10.,
    H = 15.,
    x0 = None,
    rMolecule = 0.5,
    reconstruct = False,
),
physp = nano.Params(
    Qmol = -1.,
    bulkcon = 1000.,
    dnaqsdamp = .5,
    bV = -0.1,
    rDPore = .9,
    Membraneqs = -0.0,
    ahemqs = -0.1,
    ahemuniformqs = False,
),
solverp = nano.Params(
    # TODO: add possibility of hybrid solve and PNP-only solve
    h = 1.,
    frac = 0.2,
    Nmax = 2e4,
    imax = 30,
    tol = 1e-2,
    cheapest = False,
    stokesiter = False, #True
    diffusivity_data = None,
))
defaultp = default.geop | default.physp

class Setup(solvers.Setup):
    default = default

    def __init__(self, get_geo=None, create_geo=True,
                 geop=None, physp=None, solverp=None, **params):
        if get_geo is None:
            get_geo = get_geo_default
        self.get_geo = get_geo

        self.init_params(params, geop=geop, physp=physp, solverp=solverp)
        self.init_geo(create_geo)
        self.init_phys()

    def init_geo(self, create_geo=True):
        h = self.solverp.h
        print "h", h
        self.geo = self.get_geo(h=h, **self.geop) if create_geo else None

    def init_phys(self):
        cyl = self.geop.dim == 2
        self.phys = nano.Physics("pore_mol", self.geo, cyl=cyl, **self.physp)

    def prerefine(self, visualize=False):
        return prerefine(self, visualize=visualize)


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

    # TODO:
    # solve 1D problem for side BCs
    #set_sideBCs(phys, setup.geop, setup.physp)

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
            if phys.dim==2:
                dolfin.plot(geo.boundaries, key="b", title="adapted mesh",
                            scalarbar=False)
    print "CPU time (PB): %.3g s" %(dolfin.toc(),)
    return pb

def set_D_from_data(phys, data):
    if data is not None:
        func, mesh = fields.get_functions(**data)
        D = func["D"]
        dolfin.plot(D[0], title="Dx", interactive=True)
        D = dolfin.as_matrix(np.diag([D[i] for i in range(phys.dim)]))
        phys.update(Dp=D, Dm=D)

#def solve1D(geop, physp):
#    geo = pughpore.get_geo1D(lc=.01, **geop)
#    phys = nano.Physics("pore", geo, **physp)
#    pnp = nano.solve_pde(simplepnps.SimplePNPProblem, geo, phys)
#    return geo, pnp
#
#def visualize1D(geo, pnp):
#    v, cp, cm = pnp.solutions()
#    h = geo.params["H"]
#    nano.plot1D({"potential": v}, (-h/2, h/2, 1001),
#                "x", dim=1, axlabels=("z [nm]", "potential [V]"))
#    nano.plot1D({"c+": cp, "c-":cm},  (-h/2, h/2, 1001),
#                "x", dim=1, axlabels=("z [nm]", "concentrations [mol/m^3]"))
#
#class u1D(dolfin.Expression):
#    def __init__(self, u, damping=1., **kw):
#        self.u = u
#        self.damping = damping
#        #dolfin.Expression.__init__(self)
#
#    def damp(self, scalar):
#        self.damping *= scalar
#
#    def eval(self, value, x):
#        dim = x.shape[0]
#        value[0] = self.damping*self.u(x[dim-1])
#
#def set_sideBCs(phys, geop, physp):
#    geo, pnp = solve1D(geop, physp)
#    v, cp, cm = pnp.solutions()
#    phys.v0["sideb"] = u1D(v, degree=1)
#    phys.cp0["sideb"] = u1D(cp, degree=1)
#    phys.cm0["sideb"] = u1D(cm, degree=1)

def join_dicts(list):
    "[{'F':1.0}, {'F':2.0}, ...] --> {'F':[1.0, 2.0, ...]}"
    return {key:[dic[key] for dic in list] for key in list[0]}

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

if __name__ == "__main__":
    params = nano.user_params(h=1., Nmax=1e4, dim=2, x0=[0,0,-8.5], ahemuniformqs=True)
    ddata = dict(name="Dalphahem", dim=2, Nmax=.4e5, h=1., ahemqsuniform=True, rMolecule=0.11)
    setup = Setup(**params)
    print setup.geo.params
    print setup.physp
    print setup.solverp
    setup.geo.plot_boundaries(interactive=True)
    exit()
    _, pnps = solve(setup, True)
    print get_forces(setup, pnps)

    plotter = Plotter(setup)
    v, cp, cm, u, p = pnps.solutions()
    plotter.plot_vector(u, "velocity")
    #plotter.plot(cm, "cm")
    #plotter.plot(cp, "cp")
    #plotter.plot(p, "p")
    dolfin.interactive()
