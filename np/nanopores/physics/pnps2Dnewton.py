""" Define PNP-Stokes related problems - axisymmetric version - newton method only """

from dolfin import *
from ..tools import *
from .pnps import PNPS
from .params_physical import *
from warnings import warn
from importlib import import_module

parameters["refinement_algorithm"] = "plaza_with_parent_facets"

__all__ = ["PNPSAxisymNewton","PNPSProblemAxisym"]

class PNPSAxisymNewton(PNPS):

    def __init__(self, geo, phys):
        mesh = geo.mesh
        Fmult = self.Functional_mult

        X = PNPSProblemAxisym.space(mesh)
        x = Function(X)
        v = interpolate(Constant(0.0), X.sub(0).collapse())
        c0 = interpolate(geo.pwconst("initial_ions"), X.sub(1).collapse())
        assign(x, [v, c0, c0, interpolate(Constant((0.0, 0.0)), X.sub(3).collapse()), interpolate(Constant(0.0), X.sub(4).collapse())])

        # apply BCs
        geo.BC(X.sub(0), Constant(0.), "ground").apply(x.vector())
        if phys.bV is not None:
            geo.BC(X.sub(0), Constant(phys.bV), "bV").apply(x.vector())
        geo.BC(X.sub(1), Constant(phys.bulkcon), "bulk").apply(x.vector())
        geo.BC(X.sub(2), Constant(phys.bulkcon), "bulk").apply(x.vector())
        geo.BC(X.sub(3), Constant((0.,0.)), "noslip").apply(x.vector())
        geo.BC(X.sub(4), Constant(0.), "nopressure").apply(x.vector())


        # Problem Definitions
        pnpsproblem = PNPSProblemAxisym(geo, phys, x=x)
        PNPS = IllposedNonlinearSolver(pnpsproblem)

        # Goal Functionals
        xold = pnpsproblem.solution()
        (vold, cpold, cmold, uold, pold) = xold.split()

        functionals = {}
        dim = mesh.topology().dim()
        n = FacetNormal(mesh)


        try:
            x0 = geo.parameter("x0")
        except:
            x0 = None

        if x0 is not None and "moleculeb" in geo._physical_boundary:
            # remember the following line is bad, it leads to incorrect calculation cf surface,
            # if it's not used with dS_mol(13)
            #dS_mol = geo.dS("moleculeb")
            ms_area = assemble(Expression("2*pi*x[0]")('+') * geo.dS("moleculeb"))

            Fplf = [Fmult*Expression("2*pi*x[0]")('+')*(pold)*dot(Identity(2),n('+'))[i] *geo.dS("moleculeb") for i in range(dim)]
            Fshearlf = [Fmult*Expression("2*pi*x[0]")('+')*eta*2.0*dot(sym(grad(uold))('+'),n('+'))[i] * geo.dS("moleculeb") for i in range(dim)]
            Fbarelf = [Fmult*Expression("2*pi*x[0]")('+')*qq/ms_area*grad(vold)('+')[i] * geo.dS("moleculeb") for i in range(dim)]

            F_dict = dict(
                Fp_r = Functional(Fplf[0]),
                Fp_z = Functional(Fplf[1]),
                Fshear_r = Functional(Fshearlf[0]),
                Fshear_z = Functional(Fshearlf[1]),
                Fbare_r = Functional(Fbarelf[0]),
                Fbare_z = Functional(Fbarelf[1]),
                Feff = Functional(Fplf[1] + Fshearlf[1] + Fbarelf[1]),
            )

        elif geo.parameter("name") == "P_geo":
            Fplf = [Fmult*Expression("2*pi*x[0]")('+')*(pold)*dot(Identity(2),n('+'))[i] *geo.dS("dnab") for i in range(dim)]   # pressure has opposite sign in stokes eqn
            Fshearlf = [Fmult*Expression("2*pi*x[0]")('+')*eta*2.0*dot(sym(grad(uold))('+'),n('+'))[i] * geo.dS("chargeddnab") for i in range(dim)]
            DNAqs = DNAql/geo.parameter("DNAradius")/2/pi
            f_DNAqs = Constant(DNAqs)
            Fbarelf = [Fmult*Expression("2*pi*x[0]")('+')*f_DNAqs*grad(vold)('+')[i] * geo.dS("chargeddnab") for i in range(dim)]
            Feff = Fplf[1] + Fshearlf[1] + Fbarelf[1]

            F_dict = dict(
                Fp_r = Functional(Fplf[0]),
                Fp_z = Functional(Fplf[1]),
                Fshear_r = Functional(Fshearlf[0]),
                Fshear_z = Functional(Fshearlf[1]),
                Fbare_r = Functional(Fbarelf[0]),
                Fbare_z = Functional(Fbarelf[1]),
                Feff = Functional(Feff),
            )

        else:
            F_dict = {}
        functionals.update(F_dict)


        C = geo.pwconst('diffusion_factor')
        SD = geo.pwconst('stokes_damp')

        Jm = cFarad*(C*D*grad(cmold) - C*mu*cmold*grad(vold) - SD*cmold*uold)
        Jp = cFarad*(-C*D*grad(cpold) - C*mu*cpold*grad(vold) + SD*cpold*uold)
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
            Javgtop = Functional(Expression("2*pi*x[0]")/ltop * Jz * geo.dx("poretop"))
            Jdict.update(dict(Javgtop = Javgtop))
        if "porecenter" in geo._physical_domain:
            Javgctr = Functional(Expression("2*pi*x[0]")/lcenter * Jz * geo.dx("porecenter"))
            Jdict.update(dict(Javgctr = Javgctr))
        if "porebottom" in geo._physical_domain:
            Javgbtm = Functional(Expression("2*pi*x[0]")/lbottom * Jz * geo.dx("porebottom"))
            Jdict.update(dict(Javgbtm = Javgbtm))

        functionals.update(Jdict)

        self.geo = geo
        self.phys = phys
        self.functions = {"PNPS":x}
        self.solvers = {"PNPS":PNPS}
        self.functionals = functionals

    def solve(self, refinement=False, visualize=False, save_mesh=False, print_functionals=False):
        if refinement and self.geo.mesh.num_cells() > self.maxcells:
            print 'Initial mesh has more than maximal number of cells',  \
                           ' \n  ==> no refinement \n'
            refinement = False
        tt = 0

        for i in range(self.imax):
            print '\n- Loop ' +str(i+1) + ' of max.', self.imax
            timer = Timer('Solving step '+str(i+1))

            if refinement:
                newton_iter = self.newton_solve()
                tt0 = timer.stop()
                print "Newton iterations:", newton_iter

            else:
                self.single_solve()
                tt0 = timer.stop()
                tt += tt0
                nerror = self.solvers["PNPS"].relerror()
                self.save_estimate("newton_iterations", nerror, N=i)
                self.save_estimate("newton_cputime", nerror, N=tt)
                #self.save_estimate("newton", norm(self.solvers["PNPS"].problem.u, "H10"), N=i)
                if self.solvers["PNPS"].convergence(self.tolnewton):
                    print 'linf Norm of Newton update:', \
                        norm(self.solvers["PNPS"].problem.u.vector(),'linf'), \
                        '<=', self.tolnewton ,' \n  ==> break loop \n'
                    break

            print 'Relative H1 Newton error:',\
                            self.solvers["PNPS"].relerror()

            #plot(sqrt(ind), title="sqrt(ind) "+str(i+1))
            #interactive()

            if save_mesh:
                self.save_mesh()
            if visualize:
                if visualize==True:
                    self.visualize()
                else:
                    self.visualize(visualize)
            if print_functionals:
                self.print_functionals()

            if refinement:
                (ind,err) = self.estimate()
                print "Relative error estimate (H1):",err
                self.save_estimate("h1", err)
                refined = self.refine(ind)
                if not refined:
                    tt0 = timer.stop()
                    print "Loop timing:", tt0
                    print 'Maximal number of cells reached',  \
                           ' \n  ==> no more refinement \n'
                    break

            print "Loop timing:",timer.stop()
        end


    def newton_solve(self,tol=None):
        if not tol:
            tol = self.tolnewton

        for i in range(self.imax):
            self.single_solve()
            if self.solvers["PNPS"].convergence(tol):
                break
        return i+1

    def estimate(self):
        """ simple residual indicator, estimator """
        return poisson_indicator(self.geo, self.functions["PNPS"].sub(0), cyl=True)

    def print_functionals(self):
        Jdir = self.functionals
        for Jstr in sorted(self.functionals):
            J = Jdir[Jstr]
            print ("%s: " %Jstr) + str(J.evaluate())

    def print_results(self, names=None):
        if not names:
            self.print_functionals()

    def solutions(self, string=None, deepcopy=False):
        if string:
            return self.functions[string].split(deepcopy=deepcopy)
        return self.solutions("PNPS", deepcopy)



class PNPSProblemAxisym(AdaptableNonlinearProblem):
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
        kp = "ilu",
        kparams = dict(
            maximum_iterations = 500,
            monitor_convergence = True,
            relative_tolerance = 1e-2,
            preconditioner = dict(
                ilu = dict(fill_level = 1)))
    )

    @staticmethod
    def space(mesh):
        k = PNPSProblemAxisym.k
        V = FunctionSpace(mesh, 'CG', k)
        U = VectorFunctionSpace(mesh, 'CG', k+1)
        return MixedFunctionSpace((V, V, V, U, V))

    def __init__(self, geo, phys=None, bcs=None, x=None):
        mesh = geo.mesh
        X = self.space(mesh)

        (v, cp, cm, u, p) = TrialFunctions(X)
        (vv, dp, dm, uu, q) = TestFunctions(X)

        if not x:
            x = Function(X)
            x.interpolate( Constant((0.0, phys.bulkcon, phys.bulkcon, 0.0, 0.0, 0.0)) )
            if phys.bV:
                geo.BC(X.sub(0), Constant(phys.bV), "bV").apply(x.vector())

        (vold, cpold, cmold, uold, pold) = x.split()

        dx = geo.dx()
        dx_ions = geo.dx('ions')
        r = Expression("x[0]")
        r2pi = Expression("2*pi*x[0]")
        pi = 3.14159265359

        Aperm = geo.pwconst('permittivity')
        C = geo.pwconst('diffusion_factor')
        SD = geo.pwconst('stokes_damp')

        astokes = (eta*inner(grad(u),grad(uu))*r + eta*u[0]*uu[0]/r + \
            (div(uu)*r+uu[0])*p + q*(u[0] + div(u)*r))*2*pi*dx_ions + \
            cFarad*inner((cp - cm)*grad(vold) + (cpold - cmold)*grad(v),uu)*r*2*pi*dx_ions

        Lstokes = (eta*inner(grad(uold),grad(uu))*r + eta*uold[0]*uu[0]/r + \
            (div(uu)*r+uu[0])*pold + q*(uold[0] + div(uold)*r))*2*pi*dx_ions + \
            cFarad*inner((cpold - cmold)*grad(vold),uu)*r*2*pi*dx_ions

        apoisson = inner(Aperm*grad(v),grad(vv))*r2pi*dx - cFarad*(cp - cm)*vv*r2pi*dx_ions
        aJm = inner(C*(D*grad(cm) - mu*(cm*grad(vold) + cmold*grad(v))) - SD*(cm*uold + cmold*u), grad(dm))*r2pi*dx_ions
        aJp = inner(C*(D*grad(cp) + mu*(cp*grad(vold) + cpold*grad(v))) - SD*(cp*uold + cpold*u), grad(dp))*r2pi*dx_ions
        a = astokes + apoisson + aJm + aJp

        Lpoisson = inner(Aperm*grad(vold),grad(vv))*r2pi*dx - cFarad*(cpold - cmold)*vv*r2pi*dx_ions
        LJm = inner(C*(D*grad(cmold) - mu*cmold*grad(vold)) - SD*cmold*uold, grad(dm))*r2pi*dx_ions
        LJp = inner(C*(D*grad(cpold) + mu*cpold*grad(vold)) - SD*cpold*uold, grad(dp))*r2pi*dx_ions
        Lqvol = geo.linearRHS(vv*r2pi, "volcharge")
        Lqsurf = geo.NeumannRHS(vv*r2pi, "surfcharge")
        Lq = Lqvol + Lqsurf

        L = Lstokes + Lpoisson + LJm + LJp - Lq

        if not bcs:
            try:
                bcs = [geo.BC(X.sub(0), Constant(0.0), "bV")] if phys.bV else []
                bcs = bcs + [geo.BC(X.sub(0), Constant(0.0), "ground"),
                             geo.BC(X.sub(1), Constant(0.0), "bulk"),
                             geo.BC(X.sub(2), Constant(0.0), "bulk"),
                             geo.BC(X.sub(3), Constant((0.0,0.0)), "noslip"),
                             geo.BC(X.sub(4), Constant(0.0), "nopressure"),                ]

            except:
                warning("No boundary conditions have been assigned to %s" %type(self).__name__)


        AdaptableNonlinearProblem.__init__(self, a, L, x, bcs, geo.boundaries)
