""" define solvers that record cumulative times needed for every loop """

from dolfin import *
from nanopores.tools.pdesystem import newtonsolve
from nanopores import *

def QmolEff(U, geo):
    phi, phidual = U
    phys = geo.physics
    dim = phys.dim
    coeff = Constant(1.) if dim==3 else Expression("2*pi*x[0]")
    molqv = phys.Moleculeqv
    dnaqs = phys.DNAqs
    lscale = phys.lscale
    grad = phys.grad
    q = phys.qq
    #qmol = coeff*molqv/lscale**3/q*geo.dx("molecule")
    qDNA = (1./lscale**2/q)*geo.NeumannRHS(coeff, "surfcharge")
    qmol = (1./lscale**3/q)*geo.linearRHS(coeff, "volcharge")
    Fbare = molqv * (-coeff*grad(phi)[dim-1]) *geo.dx("molecule")
    return dict(qmol=qmol, Fbare=Fbare, qDNA=qDNA)

def adaptive_pbpnps(geo, phys, cyl=False, frac=0.5, Nmax=1e4, mesh2D=None, cheapest=False,
     Felref=None, Fsref=None, Fpref=None, Fpbref=None, ratio=0.01):
    LinearPB = LinearPBAxisymGoalOriented if cyl else LinearPBGoalOriented
    PNPStokes = PNPSAxisym if cyl else PNPS
    z = phys.dim - 1
    
    bV = phys.bV
    print "biased voltage:", bV
    phys.bV = 0.
    goal = lambda v : phys.Fbare(v, z) #+ phys.Fbaresurf(v, z)
    pb = LinearPB(geo, phys, goal=goal, ref=Fpbref)
    phys.bV = bV
    pb.maxcells = Nmax
    pb.marking_fraction = frac
    if cheapest:
        pb.estimate = pb.estimate_cheap
    pb.add_functionals([QmolEff])
    refined = True
    i = 0
    
    print "Number of cells:", pb.geo.mesh.num_cells()
        
    while refined:
        i += 1
        if phys.dim == 3:
            print "\nAssessing mesh quality."
            mesh_quality(pb.geo.mesh, ratio=ratio, geo=pb.geo, plothist=False)
        print "\nSolving PB."
        # solve pb
        pb.single_solve()
        pb.print_functionals()
        
        # define and solve pnps
        '''
        if i==1:
            pnps = PNPStokes(pb.geo, phys, v0=pb.solution)
        else:
            pnps.geo = pb.geo
            mesh = pb.geo.mesh
            for name, S in pnps.solvers.items():
                print "Adapting %s." % name
                S.adapt(mesh)
            functions = tuple(pnps.functions.values())
            for S in pnps.solvers.values():
                S.replace(functions,functions)
        '''
        print "Defining PNPS with Taylor-Hood elements."
        pnps = PNPStokes(pb.geo, phys, v0=pb.solution, taylorhood=True)

        print "\nSolving PNPS."
        dofs = pnps.dofs()
        print "  Degrees of freedom: %d" % dofs
        pnps.solve()
        #newton_iter = pnps.newton_solve()
        #print "  Newton iterations:", newton_iter
        print
        #if phys.dim == 3:
        #    pnps.visualize("pore")
        fs = pnps.get_functionals()
        Fp = fs["Fp%d" %z]
        Fshear = fs["Fshear%d" %z]
        Fdrag = Fp + Fshear
        Fel = fs["Fbarevol%d" %z]
        F = Fdrag + Fel
        print "Fbare [pN]:", Fel
        print "Fdrag [pN]:", Fdrag, " = %s (Fp) + %s (Fshear)" %(Fp, Fshear)
        print "F     [pN]:", F
        if Felref is not None:
            pb.save_estimate("Fp", abs((Fp-Fpref)/Fpref), N=dofs)
            pb.save_estimate("Fel", abs((Fel-Felref)/Felref), N=dofs)
            pb.save_estimate("Fs", abs((Fshear-Fsref)/Fsref), N=dofs)
            Fref = Felref + Fsref + Fpref
            pb.save_estimate("F", abs((F-Fref)/Fref), N=dofs)
        
        print "\nAdaptive refinement."
        (ind, err) = pb.estimate()
        pb.save_estimate("Fpb est", err, N=dofs)
        refined = pb.refine(ind)
        if not refined:
            print "Maximal number of cells reached."
        else:
            print "New total number of cells:", pb.geo.mesh.num_cells()

    return pb, pnps
    
def adaptive_pb(geo, phys, cyl=False, frac=0.5, Nmax=1e4, Fpbref=None,
        ratio=.01, mesh2D=None, cheapest=False):
    LinearPB = LinearPBAxisymGoalOriented if cyl else LinearPBGoalOriented
    z = phys.dim - 1
    bV = phys.bV
    phys.bV = 0.
    goal = lambda v : phys.Fbare(v, z) #- phys.CurrentPB(v)
    pb = LinearPB(geo, phys, goal=goal, ref=Fpbref)
    phys.bV = bV
    pb.maxcells = Nmax
    pb.marking_fraction = frac
    if cheapest:
        pb.estimate = pb.estimate_cheap
    pb.add_functionals([QmolEff])
    refined = True
    i = 0
    
    print "Number of cells:", pb.geo.mesh.num_cells()
    while refined:
        i += 1
        if phys.dim == 3:
            print "\nAssessing mesh quality."
            mesh_quality(pb.geo.mesh, ratio=ratio, geo=pb.geo, plothist=False)
        print "\nSolving PB."
        # solve pb
        pb.single_solve()
        pb.print_functionals(name="Fbare")
        
        #plot(pb.geo.mesh)
        #plot(pb.geo.submesh("membrane"))
        #plot(pb.geo.submesh("pore"))
        #plot(pb.geo.submesh("dna"))
        if phys.dim == 3:
            dofs = pb.dofs()
            #plot_on_sub(pb.solution, geo, "dna", title="N=%s" %dofs)
            #geo_debug(pb.geo)
            Rz = pb.geo.params["Rz"]
            r0 = pb.geo.params["r0"]
            plot1D({"phi, N=%s" %dofs: pb.solution}, (-Rz, Rz, 101), "z", dim=3,
                origin=(r0, 0., 0.), axlabels=("z [nm]", "potential [V]"), newfig=False)
            #     origin=(0., 0., 0.), axlabels=("z [nm]", "potential [V]"), newfig=False)
        
        print "\nError estimation."
        (ind, err) = pb.estimate()
        print "\nMesh refinement."
        refined = pb.refine(ind)
        if not refined:
            print "Maximal number of cells reached."
        else:
            print "New total number of cells:", pb.geo.mesh.num_cells()
    return pb
    
def pbpnps(geo, phys, cyl=False, frac=0.5, Nmax=1e4, cheapest=False):
    LinearPB = LinearPBAxisymGoalOriented if cyl else LinearPBGoalOriented
    PNPStokes = PNPSAxisym if cyl else PNPS
    z = phys.dim - 1
    bV = phys.bV
    phys.bV = 0.
    goal = (lambda v : phys.Fbare(v, z)) if geo.parameter("x0") else (lambda v : phys.CurrentPB(v))
    pb = LinearPB(geo, phys, goal=goal)
    phys.bV = bV
    pb.maxcells = Nmax
    pb.marking_fraction = frac
    if cheapest:
        pb.estimate = pb.estimate_cheap
    refined = True
    i = 0
    
    print "Number of cells:", pb.geo.mesh.num_cells()
    while refined:
        i += 1
        print "\nSolving PB."
        pb.single_solve()
        print "\nError estimation."
        (ind, err) = pb.estimate()
        print "\nMesh refinement."
        refined = pb.refine(ind)
        if not refined:
            print "Maximal number of cells reached."
        else:
            print "New total number of cells:", pb.geo.mesh.num_cells()
            
    pnps = PNPStokes(pb.geo, phys, v0=pb.solution, taylorhood=True)
    print "\nSolving PNPS."
    dofs = pnps.dofs()
    print "  Degrees of freedom: %d" % dofs
    newton_iter = pnps.newton_solve()
    print "  Newton iterations:", newton_iter
    return pb, pnps

def newton_solve(self, tol=None, damp=None, verbose=True):
    if tol is None: tol = self.tolnewton
    if damp is None: damp = self.newtondamp
    S = self.solvers.values()[0]    
    S.newtondamp = damp
    tcum = 0.
    Uold = self.solutions(deepcopy=True)
    
    for i in range(self.imax):
        tloop = Timer("loop")
        S.solve()
        if verbose:
            print '     Relative L2 Newton error:',S.relerror()
        if S.convergence(tol):
            if verbose:
                print "     Break loop because tolerance %s was reached." %tol
            converged = True
            break
        # cumulative time 
        tcum += tloop.stop()
        # calculate the error
        U = self.solutions(deepcopy=True)
        err = sum(errornorm(u, uold, "L2", degree_rise=0) for u, uold in zip(U, Uold)) / sum(norm(u, "L2") for u in U) 
        Uold = U
        self.save_estimate("err newton i", err, N=i+1)
        self.save_estimate("err newton time", err, N=tcum)
        continue
    else:
        if verbose: print "     Did not reach tolerance %s." %tol
        converged = False
    print "     Newton iterations:",i+1
        #print '     Relative L2 Newton error:',S.relerror()
    return i+1, converged
    

def hybrid_solve(self, tol=None, damp=None):
    if tol is None: tol = self.params["tolnewton"]
    if damp is None: damp = self.params["damp"]
    I = self.params["ipicard"]
    J = self.params["inewton"]
    nverbose = self.params["nverbose"]
    verbose = self.params["verbose"]
    times = {name : 0. for name in self.solvers}
    tcum = 0.
    Uold = self.solutions(deepcopy=True)

    for i in range(1, I+1):
        #v = self.functions["npm"]
        #plot(v)
        #self.visualize()
        if verbose:
            print "\n-- Fixed-Point Loop %d of max. %d" % (i, I)
        tloop = Timer("loop")
        for name, solver in self.solvers.items():
            if verbose:
                print "    Solving %s." % name
            t = Timer(name)
            if solver.is_linear:
                solver.solve()
                times[name] += t.stop()
            else:
                j, con = newtonsolve(solver, tol, damp, J, nverbose)
                times[name] += t.stop()
                if j==1 and con:
                    print "- Break at iteration %d because Newton stopped changing." %i
                    break
        else:
            # cumulative time 
            tcum += tloop.stop()
            # calculate the error
            U = self.solutions(deepcopy=True)
            err = sum(errornorm(u, uold, "L2", degree_rise=0) for u, uold in zip(U, Uold)) / sum(norm(u, "L2") for u in U) 
            Uold = U
            self.save_estimate("err hybrid i", err, N=i)
            self.save_estimate("err hybrid time", err, N=tcum)
            continue
        break
        
    Tt = sum(times.values())
    if verbose:
        print "\n CPU Time (solve): %.2f s" % Tt
        for name in self.solvers:
            print "  -) %s: %.2f s" %(name, times[name])
            
def geo_debug(geo):
    print "Boundaries:"
    for i in geo._bou2phys:
        print "%d: %s" %(i, str(geo._bou2phys[i]))
        
    for subd in geo._physical_domain:
        submesh = geo.submesh(subd)
        geo_sub = geo_from_subdomains(submesh,
                    "nanopores.geometries.%s.subdomains" %geo.params["name"], **geo.params)
        plot(geo_sub.boundaries, title=("boundaries on %s" %subd), elevate=-3e1)
        #plot(submesh, title=("initial mesh on %s" %subd), wireframe=True, elevate=-3e1)
    interactive()
    
def mesh_quality(mesh, oldmesh=None, ratio=1e-1, geo=None, plothist=True):
    #vertex = VertexFunction("bool", mesh, False)
    dgncells = CellFunction("size_t", mesh, 0)
    
    ndeg = 0
    for c in cells(mesh):
        if c.radius_ratio() < ratio:
            dgncells[c] = 1
            ndeg += 1
    
    print "%s degenerate cells of radius ratio < %s." % (ndeg, ratio)            
    minrr = MeshQuality.radius_ratio_min_max(mesh)[0]
    print "Minimal radius ratio of mesh:", minrr
    if plothist:
        from matplotlib import pyplot
        pyplot.figure()
        exec(MeshQuality.radius_ratio_matplotlib_histogram(mesh, 200), locals())
    # plot degenerate cells
    if minrr < ratio:
        submesh = SubMesh(mesh, dgncells, 1)
        title = "degenerate N=%s" %mesh.num_cells()
        #plot(submesh, title=title)
        geo_sub = geo_from_subdomains(submesh,
                    "nanopores.geometries.%s.subdomains" %geo.params["name"], **geo.params)
        plot(geo_sub.boundaries, title="boundaries "+title)
        # find degenerate cells before snapping
        if oldmesh is not None:
            oldmesh = refine(oldmesh)
            oldcells = CellFunction("size_t", oldmesh, 0)
            oldcells.array()[:] = dgncells.array()
            plot(SubMesh(oldmesh, oldcells, 1), "old degenerate cells N=%s" %mesh.num_cells())
            
def save_Fref(pb, pnps):
    z = pnps.phys.dim - 1
    fs = pnps.get_functionals()
    Fp = fs["Fp%d" %z]
    Fshear = fs["Fshear%d" %z]
    Fdrag = Fp + Fshear
    Fel = fs["Fbarevol%d" %z]
    F = Fdrag + Fel
    Fpbref = pb.get_functionals()["goal"]
    data = dict(
        Fpref = Fp,
        Fsref = Fshear,
        Felref = Fel,
        Fpbref = Fpbref,
    )
    save_dict(data, ".", "Fref")
    
def load_Fref():
    return load_dict(".", "Fref")
    
            
            
