import numpy, dolfin
from dolfin import *
from dolfin.fem.solving import _extract_u
import petsc4py #, mpi4py
import ufl
from warnings import warn

__all__ = ["IllposedLinearSolver", "adaptform", "adaptfunction","adaptspace",
    "replace_function_in_form", "AdaptableLinearProblem",
    "IllposedNonlinearSolver",  "AdaptableNonlinearProblem",
    "AdaptableBC", "Functional", "assemble_scalar", "meshelement"]

class IllposedLinearSolver(object):
    #stabilizer constant needed for correct calculations
    stab = 1e0 #TODO ???

    def __init__(self, problem, **method):
        if not isinstance(problem, AdaptableLinearProblem):
            dolfin_error("illposed.py",
                         "initialize IllposedLinearSolver",
                         "Sorry, accepts only an AdaptableLinearProblem")
        try:
            self.method = problem.method
        except AttributeError:
            self.method = dict(reuse=False,
                               iterative=False,
                               lusolver="default",)
        if method:
            self.method.update(method)

        self.illposed = ("illposed" not in self.method) or self.method["illposed"]
        self.problem = problem

        if self.illposed:
            self.problem.a = self.stab*self.problem.a
            self.problem.L = self.stab*self.problem.L

        if self.method["reuse"]:
            self.assemble_A()

    def assemble_A(self):
        #print ("Process %s: I'm assembling a system of size %s now!" %
        #      (mpi4py.MPI.COMM_WORLD.Get_rank(), self.problem.u.function_space().dim()))
        #print "DEBUG, form:\n", self.problem.a
        A = assemble(self.problem.a, keep_diagonal=True)
        for bc in self.problem.bcs:
            bc.apply(A)
        if self.illposed:
            A.ident_zeros()

        if not self.method["iterative"]:
            self.S = LUSolver(self.method["lusolver"])
            try:
                self.S.parameters.update(self.method["luparams"])
            except KeyError:
                pass
            self.S.set_operator(A)

        else:
            if not "fieldsplit" in self.method:
                self.method["fieldsplit"] = False
            if self.method["fieldsplit"]:
                ksp = petsc4py.PETSc.KSP().create()
                ksp.setType(petsc4py.PETSc.KSP.Type.TFQMR)
                pc = ksp.getPC()
                pc.setType(petsc4py.PETSc.PC.Type.FIELDSPLIT)

                W = self.problem.u.function_space()
                is0 = petsc4py.PETSc.IS().createGeneral(W.sub(0).dofmap().dofs())
                is1 = petsc4py.PETSc.IS().createGeneral(W.sub(1).dofmap().dofs())
                pc.setFieldSplitIS(('u', is0), ('p', is1))
                pc.setFieldSplitType(0) # 0=additive
                subksps = pc.getFieldSplitSubKSP()
                subksps[0].setType("preonly")
                subksps[0].getPC().setType("hypre")
                subksps[1].setType("preonly")
                subksps[1].getPC().setType("hypre")

                A = as_backend_type(A).mat()   # export to petsc4py
                if "preconditioning_form" not in self.method:
                    ksp.setOperators(A)
                else:
                    P = assemble(self.method["preconditioning_form"], keep_diagonal=True)
                    for bc in self.problem.bcs:
                        bc.apply(P)
                    if self.illposed:
                        P.ident_zeros()
                    P = as_backend_type(P).mat()
                    ksp.setOperators(A, P)
                ksp.setFromOptions()
                self.S = ksp

            else:
                ks = (self.method["ks"] if ("ks" in self.method) else "default")
                kp = (self.method["kp"] if ("kp" in self.method) else "default")
                self.S = PETScKrylovSolver(ks, kp)
                try:
                    self.S.parameters.update(self.method["kparams"])
                except KeyError:
                    pass

                if "preconditioning_form" in self.method:
                    P = assemble(self.method["preconditioning_form"], keep_diagonal=True)
                    for bc in self.problem.bcs:
                        bc.apply(P)
                    if self.illposed:
                        P.ident_zeros()
                    self.S.set_operators(A, P)
                else:
                    self.S.set_operator(A)

    def solve(self):
        u = self.problem.u
        #plot(u.sub(0))
        if(not(self.method["reuse"])): self.assemble_A()
        b = assemble(self.problem.L)
        for bc in self.problem.bcs:
            bc.apply(b)
        #print ("Process %s: I'm solving a system of size %s now!" %
        #      (mpi4py.MPI.COMM_WORLD.Get_rank(), self.problem.u.function_space().dim()))
        if isinstance(self.S, petsc4py.PETSc.KSP):
            l = as_backend_type(b).vec()
            (x, _) = self.S.getOperators()[0].getVecs()
            print("Solving system iteratively with PETSc fieldsplit ...")
            self.S.solve(l, x)
            u.vector().set_local(x.array.astype("float_"))
        elif isinstance(self.S, PETScKrylovSolver):
            i = self.S.solve(u.vector(),b)
            # print DEBUG
            if not "kparams" in self.method or i < self.method["kparams"]["maximum_iterations"]:
                print("    PETSc Krylov solver converged in %d iterations." %i)
            else:
                print("    PETSc Krylov solver failed to converge in %d iterations." %i)
        else:
            #print "DEBUG: dim u", u.vector().array().shape
            #print "DEBUG: dim b", b.array().shape
            self.S.solve(u.vector(),b)

    def adapt(self,mesh):
        """NOTE TO SELF: if assemble isn't working after adapt,
        check if in adaptform() EVERYTHING in this form is adapted to new mesh."""
        # adapt problem to new mesh.
        # RETURNS adapted SOLUTION because it may be needed by user
        # ATTENTION: doesn't adapt functions, calling replace() may be
        # necessary as well
        self.problem.adapt(mesh)
        if(self.method["reuse"]): self.assemble_A()
        return self.problem.u

    # TODO: I have the feeling that somehow this could be done cleaner
    def damp_bcs(self, scalar):
        self.problem.damp_bcs(scalar)
        if self.method["reuse"]:
            self.assemble_A()

    def replace(self,functions,newfunctions):
        # INPUT is a tuple of the functions to be replaced by newfunctions
        # IMPORTANT: don't  supply subfunctions here, these are handled
        # automatically, i.e. if (u,v) = w.split(), only supply w.
        # TODO: reasonable interface, user input checking
        return
        for i,f in enumerate(functions):
            self.problem.a = replace_function_in_form(self.problem.a,f,newfunctions[i])
            self.problem.L = replace_function_in_form(self.problem.L,f,newfunctions[i])

    def print_method(self):
        print(self.method)
        info(self.S.parameters, True)

# END of class

class IllposedNonlinearSolver(IllposedLinearSolver):
    newtondamp = 1

    def solve(self):
        IllposedLinearSolver.solve(self)
        #plot(self.problem.u.sub(0), title="newton increment (potential)")
        self.problem.uold.vector()[:] -= self.problem.u.vector()[:]*self.newtondamp

    def adapt(self,mesh):
        IllposedLinearSolver.adapt(self,mesh)
        return self.problem.uold

    def convergence(self,tolnewton):
        if norm(self.problem.u.vector(),'linf') >= 1e12:
            warning('\n linf norm of solution: %g' %self.problem.u.vector().norm('linf'))
        return self.problem.u.vector().norm('linf') <= tolnewton

    #def relerror(self):
    #    return norm(self.problem.u,"H10")/norm(self.problem.uold,"H10")
    def relerror(self):
        norm = self.problem.uold.vector().norm('l2')
        if norm > 0.:
            return self.problem.u.vector().norm('l2')/norm
        else:
            return self.problem.u.vector().norm('l2')

class AdaptableLinearProblem(object):
    # TODO: any use for subclassing LinearVariationalProblem?
    # ATM this only adds unnessecary overhead
    # TODO: input checks, see boundaries
    def __init__(self,a,L,u,bcs=None,boundaries=None,form_compiler_parameters=None):
        self.a = a
        self.L = L
        self.u = _extract_u(u)
        self.bcs = _extract_bcs(bcs)
        self.meshid = self.u.function_space().mesh().id()
        if boundaries:
            self.boundaries = boundaries
        else:
            self.boundaries = self.bcs[0].boundaries
        #LinearVariationalProblem.__init__(self,a,L,u,bcs,form_compiler_parameters)

    def adapt(self,mesh):
        if (mesh.id() == self.meshid):
            return None
        self.meshid = mesh.id()
        V = adaptspace(self.u.function_space(),mesh)
        self.boundaries = adaptmeshfunction(self.boundaries,mesh)

        self.bcs = [bc.adapt(mesh,V,self.boundaries) for bc in self.bcs]
        self.a = adaptform(self.a,mesh)
        self.L = adaptform(self.L,mesh)
        adaptfunction(self.u, mesh, interpolate=False, assign=True)
        return V

    def damp_bcs(self, scalar): # works only for Constant PhysicalBC or Dampable bc.g
        for bc in self.bcs:
            bc.damp(scalar)

    def solution(self):
        return self.u

    def addRHS(self, L):
        self.L = self.L + L

class AdaptableNonlinearProblem(AdaptableLinearProblem):
    # TODO: compute Jacobian automatically, tolnewton etc.
    def __init__(self, a, L, uold, bcs=None, boundaries=None):
        self.uold = _extract_u(uold)
        self.bcs1 = _extract_bcs(bcs)
        bcs0 = []
        for bc in bcs:
            #print bc
            #plot(self.uold.sub(0), title="before")
            bc.apply(self.uold.vector())
            #plot(self.uold.sub(0), title="after")
            bcs0.append(bc.homogenized()) # FIXME: this severely limits bcs to PhysicalBC type

        u = Function(self.uold.function_space())
        AdaptableLinearProblem.__init__(self, a, L, u, bcs0, boundaries)

    def adapt(self, mesh):
        V = AdaptableLinearProblem.adapt(self, mesh)
        adaptfunction(self.uold, mesh, interpolate=True, assign=True)
        self.bcs1 = [bc.adapt(mesh, V, self.boundaries) for bc in self.bcs1]
        for bc in self.bcs1:
            bc.apply(self.uold.vector())

    def increment(self):
        return self.u

    def solution(self):
        return self.uold

class AdaptableBC(DirichletBC):
    # TODO: implement AdaptableBC(V,g,SubDomain)
    def __init__(self, V, g, boundaries, i, method="topological"):
        if not isinstance(boundaries,MeshFunctionSizet):
            dolfin_error("illposed.py",
                         "create AdaptableBC",
                         "Currently only implements BC defined by a MeshFunctionSizet")
        self.boundaries = boundaries
        self.i = i
        self.g = g
        DirichletBC.__init__(self, V, g, boundaries, i, method)

    def adapt(self,mesh,V=None,boundaries=None):
        # returns adapted BC
        # get adapted function space V
        # V should be provided if BC is defined on subspace
        # TODO: is there a way around assuming level 0 or 1 space component?
        if V:
            V = adaptspace(V,mesh)
            k = list(self.function_space().component())
            if k: V = V.sub(int(k[0]))
        else:
            V = adaptspace(self.function_space(),mesh)
        if not boundaries:
            boundaries = adaptmeshfunction(self.boundaries,mesh)
        self.boundaries = boundaries
        g = self.g #value()
        #print g.__class__
        if isinstance(g, Function):
            g = adaptfunction(g, mesh)
        elif isinstance(g, GenericFunction) and not isinstance(g, Constant):
            #print g.__dict__ # FIXME
            g = g.__class__()

        return AdaptableBC(V, g, self.boundaries, self.i, self.method())

    def damp(self, scalar):
        self.g.damp(scalar)

    # FIXME: this should work at least also with non-scalar Constants ...
    def damped(self, scalar):
        if isinstance(self.g, Constant):
            self.g.assign(scalar*self.g(0.))
        else:
            self.g.damp(scalar)
        V = self.function_space()
        return AdaptableBC(V, self.g, self.boundaries, self.i, self.method())


class Functional(object):
    # TODO: Functional should also be useable as linear form
    def __init__(self,form):
        self.form = form
        self.values = []

    def evaluate(self):
        try:
            e = assemble(self.form)
            #e = assemble_scalar(self.form)
        except TypeError:
            e = self.form  #maybe it's already a float
        self.values.append(e)
        return e

    def value(self):
        return self.values[-1]

    def extrapolate(self):
        return single_aitken(self.values)

    def abserror(self):
        if len(self.values) < 3:
            warn("Cannot give error estimate with fewer than 3 different evaluations")
            return 0.0
        return abs(self.value() - self.extrapolate())

    def relerror(self):
        return self.abserror()/abs(self.value())

    def adapt(self,mesh):
        self.form = adaptform(self.form,mesh)

    def replace(self,fs,newfs):
        return
        if not isinstance(fs, (list, tuple)):
            fs = (fs,)
        if not isinstance(newfs, (list, tuple)):
            newfs = (newfs,)
        for i,f in enumerate(fs):
            self.form = replace_function_in_form(self.form,f,newfs[i])

    def firstfunction(self):
        for c in self.form.coefficients():
            if isinstance(c, Function):
                f0 = c
                break
        return f0

    def __call__(self,f=None):
        if f:
            f0 = self.firstfunction()
            replace(self.form,{f0:f})
        else:
            J = self.form
        #e = assemble_scalar(J)
        e = assemble(J)
        self.values.append(e)
        return e

def meshelement(mesh):
    dim = mesh.topology().dim()
    return ["interval", "triangle", "tetrahedron"][dim-1]

def assemble_scalar(form):
    # assembles rank-0 form using the MPI communicator of the form's mesh,
    # i.e. sum values over processes only if mesh is distributed.
    # FIXME: line below does not work in dolfin >= 2016.1
    comm = form.domain().data().mpi_comm()
    z = numpy.empty(shape=(0,2), dtype="uintp")
    tl = dolfin.TensorLayout(comm, z, 0, 0, z, False)
    x = dolfin.Scalar()
    x.init(tl)
    dolfin.assemble(form, tensor=x)
    return x.get_scalar_value()

def adaptform(form, mesh):
    if dolfin.__version__ == "1.6.0":
        return adaptform_OLD(form, mesh)
    if not hasattr(form, "_compiled_form"):
        form = Form(form)
    adapted_form = adapt(form, mesh)
    adapted_form._compiled_form = form._compiled_form
    return adapted_form

def adaptform_OLD(form,mesh,adapt_coefficients=False):
    oldmesh = form.domain().data()
    if (oldmesh.id() == mesh.id()):
        return form
    # define new domain
    newdomain = mesh.ufl_domain()
    # adapt meshfunctions
    newsubdata = list(form.subdomain_data().values())[0] # assuming single domain
    for k,meshfunction in list(newsubdata.items()):
        newsubdata[k] = adaptmeshfunction(meshfunction,mesh)

    # replace domain, meshfunctions in integrals
    integrals = []
    for itg in form.integrals():
        newitg = itg.reconstruct(domain=newdomain,
                                 subdomain_data=newsubdata[itg.integral_type()])
        integrals.append(newitg)
    newform = ufl.Form(integrals)

    # replace arguments and coefficients in form
    mapping = {}
    # adapt arguments
    for argument in form.arguments():
        newargument = adaptargument(argument,mesh)
        mapping[argument] = newargument
    if (adapt_coefficients):
        # adapt coefficients
        for coeff in form.coefficients():
            adaptcoefficient(coeff,mesh) #MOD
            #newcoeff = adaptcoefficient(coeff,mesh)
            #mapping[coeff] = newcoeff

    # TODO: is there more? is there a better way to ensure everything is adapted?
    # adapt FacetNormal
    mapping[FacetNormal(oldmesh)] = FacetNormal(mesh)
    # adapt CellSize -- have to use ufl.Circumradius or replace() complains
    mapping[ufl.Circumradius(oldmesh)] = ufl.Circumradius(mesh)

    newform = replace(newform, mapping)
    return newform

def _compute_renumbering(self):
    # Include integration domains and coefficients in renumbering
    dn = self.domain_numbering()
    cn = self.coefficient_numbering()
    renumbering = {}
    renumbering.update(dn)
    renumbering.update(cn)

    # Add domains of coefficients, these may include domains not among integration domains
    k = len(dn)
    for c in cn:
        d = c.domain()
        if d is not None and d not in renumbering:
            renumbering[d] = k
            k += 1

    return renumbering

def adaptmeshfunction(meshfunction,mesh):
    #print "Meshfunction has child:", meshfunction.has_child()
    if meshfunction.has_child():
        return meshfunction.child()
    newmeshfunction = adapt(meshfunction,mesh)
    # this is to correct fenics "bug" (?) which sets some values to the highest int
    a = newmeshfunction.array()
    a[a > 1000] = 0
    return newmeshfunction

def adaptargument(argument,mesh):
    newspace = adaptspace(argument.function_space(),mesh)
    return Argument(newspace,argument.number(),part=argument.part())

def adaptspace(space, mesh):
    # only adapt if mesh is actually new
    if (space.mesh().id() == mesh.id()):
        return space
    if dolfin.__version__ == "1.6.0":
        newelement = space.ufl_element().reconstruct(domain=mesh.ufl_domain())
        return FunctionSpaceBase(mesh, newelement)
    return adapt(space, mesh)

def adaptcoefficient(coeff,mesh): #MOD
    if(isinstance(coeff, Function)):
        #coeff.assign(adaptfunction(coeff,mesh))
        return adaptfunction(coeff,mesh)
    else:
        return coeff

def adaptfunction(function, mesh, interpolate=True, assign=False):
    # important: only adapt if mesh is actually new
    if (function.function_space().mesh().id() == mesh.id()):
        return function
    newspace = adaptspace(function.function_space(), mesh)
    newfunction = Function(newspace)
    if interpolate:
        newfunction.interpolate(function)
    if assign:
        function.assign(newfunction)
    return newfunction

def replace_function_in_form(form,f,newf):
    # is supposed to replace all subfunctions of f as well
    fname = f.name()
    mapping = {}
    for c in form.coefficients():
        if (c.name() == fname):
            mapping[c] = newf
    form = replace(form,mapping)
    for i in range(f.function_space().num_sub_spaces()):
        form = replace_function_in_form(form,f.sub(i),newf.sub(i))
    return form

def _extract_bcs(bcs):
    "Extract and check argument bcs"
    if bcs is None:
        bcs = []
    elif not isinstance(bcs, (list, tuple)):
        bcs = [bcs]
    return bcs

def iterated_aitken(l):
    if not l:
        return None
    aitken = lambda l : [l[i] - (l[i+1]-l[i])**2/(l[i+2]-2*l[i+1]+l[i])
                         for i in range(len(l)-2)]
    while l:
        l0 = l
        l = aitken(l)
    return l0[-1]

def single_aitken(l):
    if len(l) < 3:
        return l[-1]
    return l[-3] - (l[-2]-l[-3])**2/(l[-1]-2*l[-2]+l[-3])
