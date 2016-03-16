""" Define PNP-Stokes related problems """

from dolfin import *
from ..tools import *
from .params_physical import *

parameters["allow_extrapolation"] = True
parameters["refinement_algorithm"] = "plaza_with_parent_facets"

__all__ = ["PNPS","PNPProblem","StokesProblem","StokesProblemEqualOrder",
           "LinearPBGoalOriented", "LinearPBProblem", "LinearPB"]

class PNPS(PDESystem):
    imax = 50
    tolnewton = 1e-3
    maxcells = 10000
    marking_fraction = 0.5
    Functional_mult = 1e12
    alwaysstokes = False

    def __init__(self, geo, phys, v0=None):
        # TODO: initialization in 3D takes more than 3 seconds, even without assembling Stokes.
        #       where is the time spent? in the imports?
        mesh = geo.mesh
        Fmult = self.Functional_mult
        StokesProblem3D = StokesProblemEqualOrder

        # set up spaces and functions
        X = PNPProblem.space(mesh)
        W = StokesProblem3D.space(mesh)

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
            #cp.vector()[:] = c0.vector()[:] - c0.vector()[:]*v.vector()[:]/phys.UT
            #cm.vector()[:] = c0.vector()[:] + c0.vector()[:]*v.vector()[:]/phys.UT
            assign(x, [v, cp, cm])
        else:
            v = interpolate(Constant(0.0), V)
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

        # scaling hack for now
        lscale = phys.lscale
        grad = phys.grad
        Fmult = Fmult
        def Cinvlscale(i):
            return Constant((1.0/lscale)**i)

        fstokes = -cFarad*(cp - cm)*grad(v)

        # Problem Definitions
        pnpproblem = PNPProblem(geo, phys, x=x, w=w)
        stokesproblem = StokesProblem3D(geo, phys, f=fstokes, w=w)

        PNP = IllposedNonlinearSolver(pnpproblem)
        Stokes = IllposedLinearSolver(stokesproblem)

        # Goal Functionals for force, current
        functionals = {}
        dim = mesh.topology().dim()
        n = FacetNormal(mesh)

        x0 = geo.params.get("x0")

        if x0 is not None:
            dS = geo.dS("moleculeb")
            dx = geo.dx("molecule")

            F_dict = {}
            for i in range(dim):
                Fp = Fmult*(-p*n[i])('-') * Cinvlscale(2)*dS
                Fshear = Fmult*(eta*2.0*dot(sym(grad(u)),-n)[i])('-') * Cinvlscale(2)*dS
                Fbare = Fmult*Constant(phys.Moleculeqs)*(-grad(v)[i])('-') * Cinvlscale(2)*dS
                Fbarevol = Fmult*Constant(phys.Moleculeqv)*(-grad(v)[i]) * Cinvlscale(3)*dx
                for F in ["Fp","Fshear","Fbare","Fbarevol"]:
                    F_dict[F+str(i)] = Functional(locals()[F])

            '''
            print "Molecule Surface Check:"
            print assemble(Constant((1.0/lscale)**2) * geo.dS("moleculeb"))
            print assemble(Constant((1.0/lscale)**2) * dS)
            print 4.*lscale**(-2)*geo.params["rMolecule"]**2*dolfin.pi

            print "Molecule Charge Check:"
            print assemble(Constant(lscale**(-2)*phys.Moleculeqs/phys.qq) * dS)
            print assemble(Constant(lscale**(-3)*phys.Moleculeqv/phys.qq) * geo.dx("molecule"))
            '''

            '''
            # remember: the following line is bad, it leads to incorrect calculation
            # of surface, if it's not used with dS_mol(13)
            #dS_mol = geo.dS("moleculeb")
            ms_area = assemble(Constant((1.0/lscale)**2) * geo.dS("moleculeb"))
            #print "Molecule surface area:",ms_area

            Fplf = [-Fmult *p*n('-')[i] *geo.dS("moleculeb") for i in range(dim)]  # pressure has opposite sign in stokes eqn
            Fshearlf = [Fmult*eta*2.0*dot(sym(grad(u)('-')),-n('-'))[i] * geo.dS("moleculeb") for i in range(dim)]
            Fbarelf = [Fmult*qq/ms_area*grad(v)('-')[i] * geo.dS("moleculeb") for i in range(dim)]
            MSlf = Fmult*Constant(1.0) * geo.dS("moleculeb")
            functionals.update({"MoleculeSurface" : Functional(MSlf)})

            F_dict = dict(
                Fp = [Functional(Fplf[i]) for i in range(dim)],
                Fshear = [Functional(Fshearlf[i]) for i in range(dim)],
                Fbare = [Functional(Fbarelf[i]) for i in range(dim)],
            ) '''
            functionals.update(F_dict)

        # currents
        C = geo.pwconst('diffusion_factor')
        SD = geo.pwconst('stokes_damp')

        Jm = cFarad*(C*(D*grad(cm) - mu*cm*grad(v)) - SD*cm*u)
        Jp = cFarad*(C*(-D*grad(cp) - mu*cp*grad(v)) + SD*cp*u)
        Jz = (Jm + Jp)[2]

        # FIXME: this is just ugly, also because geo has to have "name" in params
        """ up until now, a geo that "supports" some feature (e.g. currents) 
        had to provide the necessary information, i.e. geo.parameter("ltop") in this case.
        this is tricky for derived quantities...
        the solution is to have the derived quantity functionality in py4geo
        and always use some meta.txt as in geo_from_xml. """
        if geo.parameter("name") == "H_cyl_geo":
            ltop = (geo.parameter("l0")-geo.parameter("l1"))/2
            lctr = geo.parameter("l1")
            lbtm = ltop
        elif geo.parameter("name") == "W_3D_geo":
            l0 = geo.parameter("lsam")+geo.parameter("lau")+geo.parameter("lsin")
            lbtm = lctr = ltop = l0/3
        elif "ltop" in geo.params: # hope for the best
            ltop = geo.parameter("ltop")
            lctr = geo.parameter("lctr")
            lbtm = geo.parameter("lbtm")

        J_dict = dict(
            #Jtop = Functional(Jz('+') * Cinvlscale(2)*Fmult*geo.dS("crosstop2d")),
            #Jctrtop = Functional(Jz('+') * Cinvlscale(2)*Fmult*geo.dS("crosscentertop2d")),
            #Jctrbtm = Functional(Jz('+') * Cinvlscale(2)*Fmult*geo.dS("crosscenterbottom2d")),
            #Jbtm = Functional(Jz('+') * Cinvlscale(2)*Fmult*geo.dS("crossbottom2d")),
            Javgtop = Functional(Jz/ltop * Cinvlscale(2)*Fmult*geo.dx("poretop")),
            Javgctr = Functional(Jz/lctr * Cinvlscale(2)*Fmult*geo.dx("porecenter")),
            Javgbtm = Functional(Jz/lbtm * Cinvlscale(2)*Fmult*geo.dx("porebottom")),
            )
        functionals.update(J_dict)

        self.geo = geo
        self.phys = phys
        self.functions = {"PNP":x,"Stokes":w}
        self.solvers = {"PNP":PNP,"Stokes":Stokes}
        self.functionals = functionals

    def solve(self, refinement=False, visualize=False, save_mesh=False, print_functionals=False):
        PNPSsc = 0  # solver calls to whole PNPS System
        print "Number of cells:",self.geo.mesh.num_cells()

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
                #if self.solvers["Stokes"].problem.method["iterative"]==True and  \
                if not self.alwaysstokes and (i < 2 or  \
                    not self.solvers["PNP"].convergence(self.tolnewton*1e3) ):
                    self.solvers["PNP"].solve()
                else:
                    self.single_solve()
                    PNPSsc = PNPSsc + 1
                    
                print "Newton max error:", norm(self.solvers["PNP"].problem.u.vector(),'linf')
                #plot(self.functions["Stokes"].sub(0))

                tt0 = timer.stop()
                tt += tt0
                nerror = self.solvers["PNP"].relerror()
                self.save_estimate("fixedpoint_iterations", nerror, N=i)
                self.save_estimate("fixedpoint_cputime", nerror, N=tt)
                #self.save_estimate("fixedpoint", norm(self.solvers["PNP"].problem.u, "H10"), N=i)
                # make sure there are at least two solver call to the PNPS system
                if self.solvers["PNP"].convergence(self.tolnewton) and PNPSsc>1:
                    print 'linf Norm of Newton update:', \
                        norm(self.solvers["PNP"].problem.u.vector(),'linf'), \
                        '<=', self.tolnewton ,' \n  ==> break loop \n'
                    break

            #print 'Relative l2 Newton error:',\
            #                self.solvers["PNP"].relerror()

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
                    print "Loop timing:",tt0
                    print 'Maximal number of cells reached',  \
                           ' \n  ==> no more refinement \n'
                    break
                print "New total number of cells:",self.geo.mesh.num_cells()

            print "Loop timing:",tt0
        return i+1

    def estimate(self):
        """ simple residual indicator, estimator """
        return poisson_indicator(self.geo, self.functions["PNP"].sub(0))

    def newton_solve(self,tol=None):
        if not tol:
            tol = self.tolnewton

        for i in range(self.imax):
            #self.visualize()
            self.single_solve()
            print "Newton max error:", norm(self.solvers["PNP"].problem.u.vector(),'linf')
            #print "Newton L2 Error:", self.solvers["PNP"].relerror()
            plot(self.functions["Stokes"].sub(0))
            if self.solvers["PNP"].convergence(tol):
                break
        return i+1

    def estimate_zz(self):
        (v,cp,cm,u,p) = self.solutions()

        Aperm = self.geo.pwconst('permittivity')
        D = self.phys.D
        c0 = self.phys.bulkcon

        #fluxes, normalized to have about the same magnitude
        Jv = Aperm*grad(v)/eperm
        Jm = (D*grad(cm) - mu*cm*grad(v) - cm*u)/(D*c0)
        Jp = (D*grad(cp) + mu*cp*grad(v) - cp*u)/(D*c0)

        indv, errv = zz_indicator(v, Jv)
        indm, errm = zz_indicator(cm, Jm, self.geo.dx('fluid'))
        indp, errp = zz_indicator(cp, Jp, self.geo.dx('fluid'))

        #plot(sqrt(indv), title="sqrt(indv)")
        #plot(sqrt(indm), title="sqrt(indm)")
        #plot(sqrt(indp), title="sqrt(indp)")

        # FIXME: indicators should be multiplied, not added
        ind = Function(FunctionSpace(self.geo.mesh,"DG",0))
        ind.vector()[:] = indv.vector()[:] + indm.vector()[:] + indp.vector()[:]

        err = (errv + errm + errp)/3

        return (ind,err)

    def visualize(self, subdomain=None):
        (v,cp,cm,u,p) = self.solutions(deepcopy=True)

        mesh = self.geo.mesh
        on = ""

        if subdomain:
            on = " on " + subdomain
            mesh = self.geo.submesh(subdomain)
            for f in (v,cp,cm,u,p):
                adaptfunction(f,mesh,assign=True)

        U = Function(v.function_space())
        U.vector()[:] = (1./self.phys.UT)*v.vector()[:]
        plot(mesh,title="final mesh"+on)
        plot(U, title='el. energy [kT]'+on)
        plot(cm, title='negative ion concentration'+on)
        plot(cp, title='positive ion concentration'+on)
        plot(u, title='velocity'+on)
        plot(p, title='pressure'+on)
        interactive()

    def print_functionals(self):
        Jdir = self.functionals
        for Jstr in sorted(self.functionals):
            J = Jdir[Jstr]
            if isinstance(J,list):
                for ii in range(len(J)):
                    print ("%s[%i]: " %(Jstr,ii)) + str(J[ii].evaluate())
            else:
                print ("%s: " %Jstr) + str(J.evaluate())

    def print_results(self):
        self.print_functionals()

    def solutions(self, string=None, deepcopy=False):
        if string:
            return self.functions[string].split(deepcopy=deepcopy)
        return self.solutions("PNP", deepcopy) + self.solutions("Stokes", deepcopy)

    def rebuild(self, mesh):
        """ Assumes geometry to have geo.rebuild """
        # TODO: This is initially a lot slower than PDESystem.rebuild (why?),
        # but seems to be faster in the long run when meshes get large

        # save functional evaluations
        functionals = self.functionals
        functions = self.functions

        for f in functions:
            adaptfunction(functions[f], mesh, assign=True)

        self.geo.rebuild(mesh)
        self.__init__(self.geo)

        # TODO: is this 'hack' acceptable?
        newfunctions = {s:functions[s] for s,S in self.solvers.items()
                             if isinstance(S, IllposedNonlinearSolver)}
        oldfs = tuple(self.functions.values())
        self.functions.update(newfunctions)
        newfs = tuple(self.functions.values())

        for s,S in self.solvers.items():
            if isinstance(S, IllposedNonlinearSolver):
                # ugly
                S.problem.uold = self.functions[s]
            S.replace(oldfs,newfs)

        for Jstr,J in self.functionals.items():
            J.replace(oldfs,newfs)
            J.values = functionals[Jstr].values

    # TODO: Why is adapt not working in 3D????
    # workaround for the time being:
    #adapt = rebuild


class StokesProblem(AdaptableLinearProblem):
    k = 2
    method = dict(solvermethods.stokes)
    """
    method = dict(
        reuse = True,
        iterative = True,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = True,
            same_nonzero_pattern = True,
            reuse_factorization = True,),
        ks = "tfqmr",
        kp = "hypre_euclid",
        fieldsplit = False, #True,
        kparams = dict(
            maximum_iterations = 10000,
            monitor_convergence = False,
            # large rel.tol. together with nonzero initial guess = bad idea!!!
            relative_tolerance = 1e-8,
            # absolute tolerance must not be too large compared with newton tol
            # (but also not too low since that would be inefficient)
            absolute_tolerance = PNPS.tolnewton*1e-3,
            nonzero_initial_guess = True,
            error_on_nonconvergence = False,
            preconditioner = dict(
                report = False,
                structure = "same_nonzero_pattern",
                ilu = dict(fill_level = 1)))
    )
    """

    @staticmethod
    def space(mesh):
        k = StokesProblem.k
        U = VectorFunctionSpace(mesh, 'CG', k)
        P = FunctionSpace(mesh, 'CG', k-1)
        return U*P

    @staticmethod
    def forms(W, geo, f):
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        dx = geo.dx("fluid")

        grad = geo.physics.grad
        div = geo.physics.div
        lscale = geo.physics.lscale

        '''# scaling hack for now
        lscale = geo.parameter("nm")/nm
        def grad(u):
            return lscale*nabla_grad(u)
        def div(u):
            return lscale*transpose(nabla_div(u))
        '''
        a = (eta*inner(grad(u),grad(v)) + div(v)*p + q*div(u))*dx
        L = inner(f,v)*dx
        p = inner(grad(u), grad(v))*dx + lscale*p*q*dx
        return (a, L, p)

    def __init__(self, geo, phys, f=None, bcs=None, w=None):
        mesh = geo.mesh
        W = self.space(mesh)
        d = geo.mesh.topology().dim()
        C0 = Constant(tuple(0. for i in range(d)))
        if f is None:
            f = C0

        a, L, p = self.forms(W, geo, f)
        #self.method["preconditioning_form"] = p

        if not w:
            w = Function(W)

        if not bcs:
            try:
                bcs = [geo.BC(W.sub(0), C0, "noslip"),
                       geo.BC(W.sub(1), Constant(0.0), "nopressure")]
            except:
                warning("No boundary conditions have been assigned to %s" %type(self).__name__)

        AdaptableLinearProblem.__init__(self,a,L,w,bcs,geo.boundaries)


class PNPProblem(AdaptableNonlinearProblem):
    k = 1
    method = dict(solvermethods.bicgstab)
    """
    method = dict(
        reuse = False,
        iterative = True,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = False,
            same_nonzero_pattern = False,
            reuse_factorization = False,),
        ks = "bicgstab",
        kp = "hypre_euclid",
        kparams = dict(
            maximum_iterations = 200,
            monitor_convergence = False,
            relative_tolerance = 1e-4,
            error_on_nonconvergence = False,
            preconditioner = dict(
                ilu = dict(fill_level = 1)))
    )"""

    @staticmethod
    def space(mesh):
        V = FunctionSpace(mesh, 'CG', PNPProblem.k)
        return MixedFunctionSpace((V, V, V))

    def __init__(self, geo, phys, bcs=None, x=None, w=None):
        mesh = geo.mesh
        X = self.space(mesh)
        (v, cp, cm) = TrialFunctions(X)
        (vv, dp, dm) = TestFunctions(X)

        if not x:
            x = Function(X)
            x.interpolate(Constant((0.0, phys.bulkcon, phys.bulkcon)))
            if phys.bV:
                geo.BC(X.sub(0), Constant(phys.bV), "bV").apply(x.vector())
        if not w:
            w = Function(StokesProblem.space(mesh))

        (uold, pold) = w.split()
        (vold, cpold, cmold) = x.split()

        # scaling hack for now
        lscale = phys.lscale
        grad = phys.grad

        dx = geo.dx()
        dx_ions = geo.dx('ions')

        eps = geo.pwconst('permittivity')
        C = geo.pwconst('diffusion_factor')
        SD = geo.pwconst('stokes_damp')

        apoisson = inner(eps*grad(v),grad(vv))*dx - cFarad*(cp - cm)*vv*dx_ions
        aJm = inner(C*(D*grad(cm) - mu*(cm*grad(vold) + cmold*grad(v))) - SD*cm*uold, grad(dm))*dx_ions
        aJp = inner(C*(D*grad(cp) + mu*(cp*grad(vold) + cpold*grad(v))) - SD*cp*uold, grad(dp))*dx_ions
        a = apoisson + aJm + aJp

        Lpoisson = inner(eps*grad(vold),grad(vv))*dx - cFarad*(cpold - cmold)*vv*dx_ions
        LJm = inner(C*(D*grad(cmold) - mu*cmold*grad(vold)) - SD*cmold*uold, grad(dm))*dx_ions
        LJp = inner(C*(D*grad(cpold) + mu*cpold*grad(vold)) - SD*cpold*uold, grad(dp))*dx_ions
        Lqvol = geo.linearRHS(vv, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(vv, "surfcharge")
        Lq = Lqvol + Lqsurf

        L = Lpoisson + LJm + LJp - Lq

        '''
        Lq = Constant(0.0)*vv('+')*geo.dS('membraneb')
        try:
            x0 = geo.parameter("x0")
        except:
            x0 = None
        if x0 is not None and "moleculeb" in geo._physical_boundary:
            #ms_area = assemble(Constant((1.0/lscale)**2) * geo.dS("moleculeb"))
            #print ms_area
            #Moleculeqs = phys.Qmol/ms_area
            f_Moleculeqs = Constant(phys.Moleculeqs)
            LqMolecule = f_Moleculeqs*vv('+')*geo.dS('moleculeb')
            Lq = Lq + LqMolecule

        if "chargeddnab" in geo._physical_boundary:
            #if "DNAqs" in phys_params:
            #    DNAqs = phys_params["DNAqs"]
            #else:
            #    DNAQ = geo.parameter("ncbp")*bpq
            #    cDNA_area = assemble(Constant((1.0/lscale)**2) * geo.dS("chargeddnab"))
            #    DNAqs = DNAQ/cDNA_area*self.dnaqsdamp
            f_DNAqs = Constant(phys.DNAqs)
            LqDNA = f_DNAqs*vv('+')*geo.dS('chargeddnab')
            Lq = Lq + LqDNA

        if "chargedsinb" in geo._physical_boundary:
            f_SiNqs = Constant(phys_params["SiNqs"])
            LqSiN = f_SiNqs*vv('+')*geo.dS('chargedsinb')
            Lq = Lq + LqSiN
        if "chargedsamb" in geo._physical_boundary:
            f_SAMqs = Constant(phys_params["SAMqs"])
            LqSAM = f_SAMqs*vv('+')*geo.dS('chargedsamb')
            Lq = Lq + LqSAM

        if "membraneb" in geo._physical_boundary:
            f_Membraneqs = Constant(phys.Membraneqs)
            LqMembrane = f_Membraneqs*vv('+')*geo.dS('membraneb')
            Lq = Lq + LqMembrane

        L = Lpoisson + LJm + LJp - lscale*Lq
        '''
        if not bcs:
            try:
                bcs = [geo.BC(X.sub(0), Constant(0.0), "bV")] if phys.bV else []
                bcs += [geo.BC(X.sub(0), Constant(0.0), "ground"),
                        geo.BC(X.sub(1), Constant(0.0), "bulk"),
                        geo.BC(X.sub(2), Constant(0.0), "bulk")]
            except:
                warning("No boundary conditions have been assigned to %s" %type(self).__name__)

        AdaptableNonlinearProblem.__init__(self, a, L, x, bcs, geo.boundaries)


class StokesProblemEqualOrder(StokesProblem):
    k = 1
    # stabilization parameter
    beta = 0.01

    @staticmethod
    def space(mesh):
        k = StokesProblemEqualOrder.k
        U = VectorFunctionSpace(mesh, 'CG', k)
        P = FunctionSpace(mesh, 'CG', k)
        return U*P

    @staticmethod
    def forms(W, geo, f):
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        dx = geo.dx("fluid")
        mesh = geo.mesh

        grad = geo.physics.grad
        div = geo.physics.div
        lscale = geo.physics.lscale

        h = CellSize(mesh)/lscale
        delta = StokesProblemEqualOrder.beta*h**2

        # added stabilization term
        a = (2*eta*inner(sym(grad(u)), sym(grad(v))) + div(v)*p + q*div(u))*dx \
             - delta*inner(grad(p),grad(q))*dx
        L = inner(f,v - delta*grad(q))*dx
        p = 2*inner(sym(grad(u)), sym(grad(v)))*dx + lscale*inner(p, q)*dx #- delta*inner(grad(p),grad(q))*dx
        return (a, L, p)

from .poisson import PoissonProblem
class LinearPBProblem(PoissonProblem):

    @staticmethod
    def forms(V, geo, f):
        u = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
        dx0 = geo.dx("ions")
        c0 = geo.physics.bulkcon
        UT = geo.physics.UT
        grad = geo.physics.grad
        lscale = geo.physics.lscale

        eps = geo.pwconst('permittivity')
        a = eps*inner(grad(u), grad(v))*dx  + Constant(cFarad*2*c0/UT)*u*v*dx0
        Lqvol = geo.linearRHS(v, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(v, "surfcharge")
        Lq = Lqvol + Lqsurf
        L = f*v*dx + Lq

        return (a, L)

class LinearPBGoalOriented(GoalAdaptivePDE):
    def __init__(self, geo, phys, goal=None, ref=None):
        if goal is None and geo.params["x0"]:
            goal = lambda v : phys.Fbare(v, 2)
        self.ref = ref # reference value for functional
        GoalAdaptivePDE.__init__(self, geo, phys, LinearPBProblem, goal)

    def estimate(self):
        u = self.functions["primal"]
        z = self.functions["dual"]
        ind, err, rep, errc, gl, glx = pb_indicator_GO(self.geo, self.phys, u, z)
        self.save_estimate("err", err)
        self.save_estimate("rep", rep)
        self.save_estimate("goal", gl)
        self.save_estimate("goal ex", glx)
        return ind, rep

    def estimate_cheap(self):
        u = self.functions["primal"]
        z = self.functions["dual"]
        ind, err, gl = pb_indicator_GO_cheap(self.geo, self.phys, u, z)
        self.save_estimate("err", err)
        self.save_estimate("goal", gl)
        return ind, err

    def print_functionals(self, name="goal"):
        PDESystem.print_functionals(self)
        J = self.functionals[name]
        Jval = J.value()
        if self.ref is not None:
            ref = self.ref
            err = abs((Jval-ref)/ref)
            self.save_estimate("err ref", err)
            self.save_estimate("goal ref", ref)

class LinearPB(LinearPDE):
    def __init__(self, geo, phys):
        LinearPDE.__init__(self, geo, LinearPBProblem, phys=phys)

    def estimate(self):
        c0 = self.geo.physics.bulkcon
        u = self.functions.values()[0]
        chi = self.geo.pwconst("ions", value={"ions":1.,"solid":0.})
        UT = self.geo.physics.UT
        f = Constant(-cFarad*2*c0/UT)*u*chi
        ind,err = poisson_indicator(self.geo, u, f=f)
        self.save_estimate("err", err)
        return ind, err
