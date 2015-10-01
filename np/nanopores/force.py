#!/usr/bin/env python

''' postprocess data file from functional calculations '''

def post_process(visualize=False, datafile=None):
    # RETURNS F = Fp + Fshear + Fbare as a function of x=[r,z]
    from nanopores import DATADIR
    from dolfin import Mesh, CellFunction, FunctionSpace, Function, cells, \
                       plot, interactive, VectorFunctionSpace, assign, parameters
    import sys
        
    parameters["allow_extrapolation"] = True

    # which functionals you want to post-process: (subset of pnps.functionals and fkeys in simulation_script)
    fkeys = {"Javgtop","Javgctr","Javgbtm","Fp0","Fp1","Fp2","Fshear0",
        "Fshear1","Fshear2","Fbare0","Fbare1","Fbare2","Fbarevol0","Fbarevol1","Fbarevol2"}
    
    # --- get filenames and meshes
    
    # get git repo
    with open("%s/sim/gitrepo"%(DATADIR,), "r") as f:
        gitrepo = f.next().strip()
    
    # get data filename
    if datafile is None:
        with open("%s/sim/datafile"%(DATADIR,), "r") as f:
            datafile = f.next().strip()
    
    folder = "np-sim/simulations"
    
    # get meshes and simulation data    
    FILENAME = "/".join([gitrepo, folder, datafile])
    mesh = Mesh(".".join(FILENAME.split(".")[:-1]) + "_mesh_i.xml")
    mesh_o = Mesh(".".join(FILENAME.split(".")[:-1]) + "_mesh_o.xml")
    
    # insert gitrepo to path and import protocol TODO?
    sys.path.insert(0, gitrepo+"/np-sim")
    from protocol import Data
    
    # --- access data file
    data = Data(FILENAME)
    
    # check whether all vertices have status 1
    calculation_complete = all(data.data["status"] == 1)
    if not calculation_complete:
        if visualize:
            # create indicator function of completed cells
            status = iter(data.data["status"])
            calc = CellFunction("size_t", mesh, 0)
            for c in cells(mesh):
                if status.next() == 1:
                    calc[c] = 1
            plot(calc, title="Completed cells")
            interactive()
        raise SystemExit("""Error: Calculation in %s is not complete.
Consider changing the simulation file with np-sim init.""" %datafile)

    # get functionals data
    functional_iter = {}
    for k in fkeys:
        functional_iter[k] = iter(data.data[k])

    # --- assign functionals to Functions
    
    # initialize
    functionals = {}
    D = FunctionSpace(mesh, "DG", 0)
    Ddofmap = D.dofmap()
    
    for k in fkeys:
        functionals[k] = Function(D)
    current = Function(D)
        
    # read from data to DG Functions
    for c in cells(mesh):
        i = Ddofmap.cell_dofs(c.index())
        for k in functionals:
            functionals[k].vector()[i] = functional_iter[k].next()
        if c.midpoint().y() < 0:
            current.vector()[i] = functionals["Javgtop"].vector()[i]
        else:
            current.vector()[i] = functionals["Javgbtm"].vector()[i]
    functionals["Jz"] = current
            
    # interpolate to CG Functions
    V = FunctionSpace(mesh_o, "CG", 1)
    for k in functionals:
        f = Function(V)
        f.interpolate(functionals[k])
        f.set_allow_extrapolation(False)
        functionals[k] = f
        
    parameters["allow_extrapolation"] = False
        
    for i in range(3):
        f = Function(V)
        f.assign(sum(functionals["%s%d"%(s,i)] for s in ["Fp", "Fshear", "Fbarevol"]))
        functionals["F%d"%i] = f

    # assign forces to vector valued functions     
    vforces = {}
    VV = VectorFunctionSpace(mesh_o, "CG", 1)
    
    for k in {"Fp", "Fshear", "Fbare", "Fbarevol", "F"}:
        functionals[k] = [functionals["%s%d"%(k,j)] for j in range(3)]
        vforces[k] = Function(VV)
        assign(vforces[k], [functionals[k+"0"], functionals[k+"2"]])
        
    fmaps = {}
    
    # Functionals can return some default value when argument is out of range:
    def make_functional_map(k, default=[0.,0.,0.]):
        #return lambda x: map(lambda j : functionals[k][j](x), [0,1,2])
        def thefunction(x):
            try:
                return map(lambda j : functionals[k][j](x), [0,1,2])
            except RuntimeError:
                print "Warning: %s lies outside of domain!" %str(x)
                return default
        return thefunction
        
    def make_functional_scalar(k, default=0.):
        #return lambda x: map(lambda j : functionals[k][j](x), [0,1,2])
        def thefunction(x):
            try:
                return functionals[k](x)
            except RuntimeError:
                return default
        return thefunction
            
    fmaps = {k:make_functional_map(k) for k in {"Fp", "Fshear", "Fbare", "F"}}
    
    #F = lambda x: map(lambda j : functionals["F"][j](x), [0,1,2])
    Jz = make_functional_scalar("Jz")
    F = make_functional_map("F")
    
    if visualize:
        x = cells(mesh).next().midpoint()
        s = "F at %s: %s"
        print s % ([x[0],x[1]], F(x))
        print s % ([0.0,0.0], F([0.0, 0.0]))
        print s % ([0.0,1e-7], F([0.0, 1e-7]))
        x = [4.910359738251182e-10, 5.132699214904982e-09]
        print s % (x, F(x))
        print "Jz at %s: %s" % ([0.0,1e-7], Jz([0.0, 1e-7]))
      
        for k in vforces:
            plot(vforces[k], title=k)
        plot(functionals["Jz"], title="Jz")
        plot(functionals["F2"], title="F2")
        interactive()
       
    return F,Jz,fmaps,functionals

if __name__ == "__main__":
    F, Jz, fmaps, functionals = post_process(visualize=True)
else:
    F, Jz, fmaps, functionals = post_process()
