""" Define PNP-Stokes related problems - axisymmetric version """

from dolfin import *
from ..tools import *
from .pnps import PNPS
from .params_physical import *
from warnings import warn
from importlib import import_module
from .simplepnps import SimpleStokesProblem

parameters["refinement_algorithm"] = "plaza_with_parent_facets"

__all__ = ["PNPSAxisym","PNPProblemAxisym","StokesProblemAxisym","LinearPBAxisym",
           "LinearPBAxisymGoalOriented", "StokesProblemAxisymEqualOrder"]

class PNPSAxisym(PNPS):

    def __init__(self, geo, phys, v0=None):
        mesh = geo.mesh
        Fmult = self.Functional_mult
        StokesProblem2D = StokesProblemAxisym #StokesProblemAxisymEqualOrder

        # set up spaces and functions
        X = PNPProblemAxisym.space(mesh)
        W = StokesProblem2D.space(mesh)

        x = Function(X)
        w = Function(W)

        V = X.sub(0).collapse()
        c0 = interpolate(geo.pwconst("initial_ions"), V)

        # optional initial guess for potential from PB eq
        if v0 is not None:
            v = interpolate(v0, V)
            import numpy
            cp = Function(V)
            cm = Function(V)
            cp.vector()[:] = c0.vector()[:]*numpy.exp(-v.vector()[:]/phys.UT)
            cm.vector()[:] = c0.vector()[:]*numpy.exp(v.vector()[:]/phys.UT)
            #cp = project(c0*exp(-v0/phys.UT), V)
            #cm = project(c0*exp(v0/phys.UT), V)
            assign(x, [v, cp, cm])
        else:
            v = interpolate(Constant(0.0),V)
            assign(x, [v, c0, c0])

        # apply BCs
        geo.BC(X.sub(0), Constant(0.), "ground").apply(x.vector())
        if phys.bV is not None:
            geo.BC(X.sub(0), Constant(phys.bV), "bV").apply(x.vector())
        geo.BC(X.sub(1), Constant(phys.bulkcon), "bulk").apply(x.vector())
        geo.BC(X.sub(2), Constant(phys.bulkcon), "bulk").apply(x.vector())

        # apply concentration bias, e.g. upperpcon = 1, lowerpcon = -1
        for side in ["upper","lower"]:
            for i, sign in enumerate(["p","m"]):
                bias = getattr(phys, side + sign + "bias")
                if bias is not None:
                    con = Constant(phys.bulkcon + bias)
                    geo.BC(X.sub(i+1), con, side+"b").apply(x.vector())

        (v, cp, cm) = x.split()
        (u, p) = w.split()

        grad = phys.grad
        fstokes = -cFarad*(cp - cm)*grad(v)

        # Problem Definitions
        pnpproblem = PNPProblemAxisym(geo, phys, x=x, w=w)
        stokesproblem = StokesProblem2D(geo, phys, f=fstokes, w=w)

        PNP = IllposedNonlinearSolver(pnpproblem)
        Stokes = IllposedLinearSolver(stokesproblem)

        # Goal Functionals
        # FIXME: this is a mess -- functionals need to be outsourced
        functionals = {}
        F_dict = {}
        dim = mesh.topology().dim()
        n = FacetNormal(mesh)
        r2pi = Expression("2*pi*x[0]")
        def tang(x):
            return x - inner(x,n)*n

        try:
            x0 = geo.parameter("x0")
        except:
            x0 = None

        if x0 is not None and "moleculeb" in geo._physical_boundary:
            dS = geo.dS("moleculeb")
            dx = geo.dx("molecule")
            rho = phys.Moleculeqs

            for i in range(dim):
                Fp = Fmult*(-r2pi*p*n[i])('-') *dS
                Fshear = Fmult*(r2pi*eta*2.0*dot(sym(grad(u)),-n)[i])('-') * dS
                Fbare = Constant(Fmult*rho)*(-r2pi*grad(v)[i])('-') * dS
                Fbarevol = Constant(Fmult*phys.Moleculeqv)*(-r2pi*grad(v)[i]) * dx

                for F in ["Fp","Fshear","Fbare", "Fbarevol"]:
                    F_dict[F+str(i)] = Functional(locals()[F])


        if geo.parameter("name") == "P_geo":
            Fplf = [Fmult*r2pi('+')*(p)*dot(Identity(2),n('+'))[i] *geo.dS("dnab") for i in range(dim)]   # pressure has opposite sign in stokes eqn
            Fshearlf = [Fmult*r2pi('+')*eta*2.0*dot(sym(grad(u))('+'),n('+'))[i] * geo.dS("chargeddnab") for i in range(dim)]

            f_DNAqs = Constant(phys.DNAqs)
            Fbarelf = [Fmult*r2pi('+')*f_DNAqs*grad(v)('+')[i] * geo.dS("chargeddnab") for i in range(dim)]
            Feff = Fplf[1] + Fshearlf[1] + Fbarelf[1]

            F_dict = dict(
                Fp0 = Functional(Fplf[0]),
                Fp1 = Functional(Fplf[1]),
                Fshear0 = Functional(Fshearlf[0]),
                Fshear1 = Functional(Fshearlf[1]),
                Fbare0 = Functional(Fbarelf[0]),
                Fbare1 = Functional(Fbarelf[1]),
                Feff = Functional(Feff),
            )

        functionals.update(F_dict)

        # currents
        C = geo.pwconst('diffusion_factor')
        SD = geo.pwconst('stokes_damp')
        Jm = cFarad*(C*(D*grad(cm) - mu*cm*grad(v)) - SD*cm*u)
        Jp = cFarad*(C*(-D*grad(cp) - mu*cp*grad(v)) + SD*cp*u)
        Jz = Fmult*(Jp + Jm)[1]

        if geo.parameter("name") == "P_geo" or geo.parameter("name") == "H_geo":
            ltop = -(geo.parameter("l1") - geo.parameter("l0"))/2
            lcenter = geo.parameter("l1")
            lbottom = ltop
        elif geo.parameter("name") == "W_2D_geo":
            ltop = (geo.parameter("lsin")+geo.parameter("lau")+geo.parameter("lsam"))/3
            lcenter = lbottom = ltop

        Jdict = {}
        if "poretop" in geo._physical_domain:
            Javgtop = Functional(r2pi/ltop * Jz * geo.dx("poretop"))
            Jdict.update(dict(Javgtop = Javgtop))
        if "porecenter" in geo._physical_domain:
            Javgctr = Functional(r2pi/lcenter * Jz * geo.dx("porecenter"))
            Jdict.update(dict(Javgctr = Javgctr))
        if "porebottom" in geo._physical_domain:
            Javgbtm = Functional(r2pi/lbottom * Jz * geo.dx("porebottom"))
            Jdict.update(dict(Javgbtm = Javgbtm))

        functionals.update(Jdict)

        self.geo = geo
        self.phys = phys
        self.functions = {"PNP":x,"Stokes":w}
        self.solvers = {"PNP":PNP,"Stokes":Stokes}
        self.functionals = functionals

    # FIXME: currently does not work because of some fenics segfault-type bug
    @staticmethod
    def prerefine(geo, phys, N=None):
        ''' cheap method to determine adapted mesh by solving only poisson-like equation.
        obviously, this does not yield any error estimate on functionals of relevance.
        needs only Geometry and Physics as input '''
        poisson = LinearPBAxisym(geo, phys)
        poisson.maxcells = PNPSAxisym.maxcells if N is None else N

        poisson.solve(refinement=True)
        #plot(poisson.solutions()[0])
        #plot(poisson.geo.mesh)
        #plot(poisson.geo.subdomains)
        #plot(poisson.geo.boundaries)
        #interactive()
        return poisson

    def estimate(self):
        """ simple residual indicator, estimator """
        return poisson_indicator(self.geo, self.functions["PNP"].sub(0), cyl=True)

    def print_functionals(self):
        Jdir = self.functionals
        for Jstr in sorted(self.functionals):
            J = Jdir[Jstr]
            print ("%s: " %Jstr) + str(J.evaluate())

    def print_results(self, names=None):
        if not names:
            self.print_functionals()

class StokesProblemAxisym(AdaptableLinearProblem):
    k = 2
    method = dict(
        reuse = True,
        iterative = False,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = True,
            same_nonzero_pattern = True,
            reuse_factorization = True,),
        ks = "minres",
        kp = "ilu",
        kparams = dict(
            maximum_iterations = 500,
            monitor_convergence = False,
            relative_tolerance = 1e-2,
            nonzero_initial_guess = True,
            preconditioner = dict(
                report = False,
                structure = "same_nonzero_pattern",
                ilu = dict(fill_level = 1)))
    )

    @staticmethod
    def space(mesh):
        k = StokesProblemAxisym.k
        U = VectorFunctionSpace(mesh, 'CG', k)
        P = FunctionSpace(mesh, 'CG', k-1)
        return U*P

    @staticmethod
    def forms(W, geo, f):
        mesh = geo.mesh

        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)

        dx = geo.dx("fluid")
        r = Expression("x[0]")
        pi = 3.14159265359
        grad = geo.physics.grad
        div = geo.physics.div

        # conservative formulation for correct BC
        a = (2*eta*inner(sym(grad(u)),sym(grad(v)))*r + 2*eta*u[0]*v[0]/r + \
            (div(v)*r+v[0])*p + q*(u[0] + div(u)*r))*2*pi*dx

        L = 2*pi*inner(f,v)*r*dx if f else inner(Constant((0.0,0.0)),v)*dx
        return (a,L)

    def __init__(self, geo, phys=None, f=None, bcs=None, w=None):
        mesh = geo.mesh
        W = self.space(mesh)
        a, L = self.forms(W, geo, f)

        if not w:
            w = Function(W)

        if not bcs:
            try:
                if "inflow" in geo._physical_boundary:
                    bcs = [geo.BC(W.sub(0), Constant((0.0,0.0)), "noslip"),
                           geo.BC(W.sub(0), Constant(phys.inflow), "inflow"),
                           geo.BC(W.sub(1), Constant(0.0), "nopressure")]
                else:
                    bcs = [geo.BC(W.sub(0), Constant((0.0,0.0)), "noslip"),
                           geo.BC(W.sub(1), Constant(0.0), "nopressure")]
            except:
                warning("No boundary conditions have been assigned to %s" %type(self).__name__)

        AdaptableLinearProblem.__init__(self, a, L, w, bcs, geo.boundaries)

class PNPProblemAxisym(AdaptableNonlinearProblem):
    # TODO: automatically differentiate nonlinear form
    k = 1
    method = dict(
        reuse = False,
        iterative = False,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = False,
            same_nonzero_pattern = False,
            reuse_factorization = False,),
        ks = "bicgstab",
        kp = "hypre_euclid",
        kparams = dict(
            maximum_iterations = 1000,
            monitor_convergence = True,
            relative_tolerance = 1e-2,
            preconditioner = dict(
                ilu = dict(fill_level = 1)))
    )

    @staticmethod
    def space(mesh):
        V = FunctionSpace(mesh, 'CG', PNPProblemAxisym.k)
        return MixedFunctionSpace((V, V, V))

    def __init__(self, geo, phys=None, bcs=None, x=None, w=None):
        mesh = geo.mesh
        V = FunctionSpace(mesh, 'CG', self.k)
        X = self.space(mesh)

        (v, cp, cm) = TrialFunctions(X)
        (vv, dp, dm) = TestFunctions(X)

        if not x:
            x = Function(X)
            x.interpolate(Constant((0.0, phys.bulkcon, phys.bulkcon)))
            if phys.bV:
                geo.BC(X.sub(0), Constant(phys.bV), "bV").apply(x.vector())
        if not w:
            w = Function(StokesProblemAxisym.space(mesh))

        (uold, pold) = w.split()
        (vold, cpold, cmold) = x.split()

        dx = geo.dx()
        dx_ions = geo.dx('ions')
        r2pi = Expression("2*pi*x[0]")
        grad = phys.grad
        lscale = Constant(phys.lscale)

        Aperm = geo.pwconst('permittivity')
        C = geo.pwconst('diffusion_factor')
        SD = geo.pwconst('stokes_damp')

        apoisson = inner(Aperm*grad(v),grad(vv))*r2pi*dx - cFarad*(cp - cm)*vv*r2pi*dx_ions
        aJm = inner(C*(D*grad(cm) - mu*(cm*grad(vold) + cmold*grad(v))) - SD*cm*uold, grad(dm))*r2pi*dx_ions
        aJp = inner(C*(D*grad(cp) + mu*(cp*grad(vold) + cpold*grad(v))) - SD*cp*uold, grad(dp))*r2pi*dx_ions
        a = apoisson + aJm + aJp

        Lpoisson = inner(Aperm*grad(vold),grad(vv))*r2pi*dx - cFarad*(cpold - cmold)*vv*r2pi*dx_ions
        LJm = inner(C*(D*grad(cmold) - mu*cmold*grad(vold)) - SD*cmold*uold, grad(dm))*r2pi*dx_ions
        LJp = inner(C*(D*grad(cpold) + mu*cpold*grad(vold)) - SD*cpold*uold, grad(dp))*r2pi*dx_ions
        Lqvol = geo.linearRHS(vv*r2pi, "volcharge")
        Lqsurf = geo.NeumannRHS(lscale*vv*r2pi, "surfcharge")
        Lq = Lqvol + Lqsurf

        L = Lpoisson + LJm + LJp - Lq
        
        if not bcs:
            try:
                bcs = [geo.BC(X.sub(0), Constant(phys.bV), "bV")] if phys.bV else []
                bcs += [geo.BC(X.sub(0), Constant(0.0), "ground"),
                        geo.BC(X.sub(1), Constant(phys.bulkcon), "bulk"),
                        geo.BC(X.sub(2), Constant(phys.bulkcon), "bulk")]
            except:
                warning("No boundary conditions have been assigned to %s" %type(self).__name__)
        """
        if not bcs:
            try:
                bcs = [geo.BC(X.sub(0), Constant(0.0), "bV")] if phys.bV else []
                bcs += [geo.BC(X.sub(0), Constant(0.0), "ground"),
                        geo.BC(X.sub(1), Constant(0.0), "bulk"),
                        geo.BC(X.sub(2), Constant(0.0), "bulk")]
            except:
                warning("No boundary conditions have been assigned to %s" %type(self).__name__)
        """
        AdaptableNonlinearProblem.__init__(self, a, L, x, bcs, geo.boundaries)


from .poisson import PoissonProblem
class PoissonProblemAxisym(PoissonProblem):
    @staticmethod
    def forms(V, geo, f):
        u = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
        r2pi = Expression("2*pi*x[0]")
        grad = geo.physics.grad
        lscale = Constant(geo.physics.lscale)

        eps = geo.pwconst('permittivity')
        a = inner(eps*grad(u), grad(v))*r2pi*dx
        Lqvol = geo.linearRHS(v*r2pi, "volcharge")
        Lqsurf = geo.NeumannRHS(lscale*v*r2pi, "surfcharge")
        Lq = Lqvol + Lqsurf
        L = f*v*r2pi*dx + Lq
        return (a, L)

class PoissonAxisym(LinearPDE):
    def __init__(self, geo, phys):
        LinearPDE.__init__(self, geo, PoissonProblemAxisym, phys=phys)
    def estimate(self):
        return poisson_indicator(self.geo, self.functions.values()[0], cyl=True)

class LinearPBProblemAxisym(PoissonProblem):

    @staticmethod
    def forms(V, geo, f):
        u = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
        dx0 = geo.dx("ions")
        r2pi = Expression("2*pi*x[0]")
        c0 = geo.physics.bulkcon
        UT = geo.physics.UT
        grad = geo.physics.grad
        lscale = Constant(geo.physics.lscale)
        print "DEBUG:", geo.physics.lscale

        eps = geo.pwconst('permittivity')
        a = eps*inner(grad(u), grad(v))*r2pi*dx  + Constant(cFarad*2*c0/UT)*u*v*r2pi*dx0
        Lqvol = geo.linearRHS(v*r2pi, "volcharge")
        Lqsurf = geo.NeumannRHS(lscale*v*r2pi, "surfcharge")
        Lq = Lqvol + Lqsurf
        L = f*v*r2pi*dx + Lq

        return (a, L)

class LinearPBAxisym(LinearPDE):
    def __init__(self, geo, phys):
        LinearPDE.__init__(self, geo, LinearPBProblemAxisym, phys=phys)

    def estimate(self):
        u = self.functions.values()[0]
        ind,err = pb_indicator(self.geo, self.geo.physics, u, cyl=True)
        self.save_estimate("err", err)
        return ind, err

class LinearPBAxisymGoalOriented(GoalAdaptivePDE):
    def __init__(self, geo, phys, goal=None):
        if goal is None:
            goal = lambda v : phys.Fbare(v, 1)
        GoalAdaptivePDE.__init__(self, geo, phys, LinearPBProblemAxisym, goal)

    def estimate(self):
        u = self.functions["primal"]
        z = self.functions["dual"]
        ind, err = pb_indicator_GO(self.geo, self.phys, u, z, cyl=True)
        self.save_estimate("err", err)
        return ind, err

    def estimate0(self):
        # primal and dual indicators
        (indp, errp) = pb_indicator(self.geo, self.phys, self.functions["primal"], cyl=True)
        (indd, errd) = pb_indicator(self.geo, self.phys, self.functions["dual"], cyl=True)

        # goal-oriented indicator
        W = FunctionSpace(self.geo.mesh, "DG", 0)
        ind = Function(W)
        import numpy
        ind.vector()[:] = numpy.sqrt(indp.vector()[:]*indd.vector()[:])
        err = sqrt(errp*errd)

        # plot indicators
        '''
        plotind = Function(W)
        plotind.vector()[:] = numpy.log(indp.vector()[:])
        plot(plotind, title="primal indicators")
        plotind.vector()[:] = numpy.log(indd.vector()[:])
        plot(plotind, title="dual indicators")
        plotind.vector()[:] = numpy.log(ind.vector()[:])
        plot(plotind, title="goal-oriented indicators")
        interactive()
        '''
        return ind, err

    def print_functionals(self):
        J = self.functionals["goal"]
        # self.save_estimate("Fbare", 1e12*self.functionals["goal"].value()) # FIXME: error: list index out of range
        print "Goal (*1e12):", J()*1e12


class StokesProblemAxisymEqualOrder(StokesProblemAxisym):
    k = 1
    beta = 0.1

    @staticmethod
    def space(mesh):
        k = StokesProblemAxisymEqualOrder.k
        U = VectorFunctionSpace(mesh, 'CG', k)
        P = FunctionSpace(mesh, 'CG', k)
        return U*P

    @staticmethod
    def forms(W, geo, f):
        mesh = geo.mesh

        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)

        dx = geo.dx("fluid")
        grad = geo.physics.grad
        div = geo.physics.div
        lscale = geo.physics.lscale
        
        r = Expression("x[0]")
        pi = 3.14159265359
        h = CellSize(mesh)
        
        beta = StokesProblemAxisymEqualOrder.beta
        print "DEBUG: beta = ", beta
        print "DEBUG: lscale = ", lscale
        delta = Constant(beta/lscale**2)*h**2

        # conservative formulation for correct BC, with added stabilization term
        a = (2*eta*inner(sym(grad(u)),sym(grad(v)))*r + 2*eta*u[0]*v[0]/r + \
            (div(v)*r+v[0])*p + q*(u[0] + div(u)*r))*2*pi*dx - \
            delta*inner(grad(p),grad(q))*r*2*pi*dx

        if f is None:
            f = Constant((0.0,0.0))

        L = 2*pi*inner(f,v - delta*grad(q))*r*dx
        return (a, L)
