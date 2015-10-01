#!/usr/bin/env python

''' postprocess data file from force calculations '''

def post_process(visualize=False):
    # RETURNS F = Fp + Fshear + Fbare as a function of x=[r,z]
    from nanopores import DATADIR
    from dolfin import Mesh, CellFunction, FunctionSpace, Function, cells, \
                       plot, interactive, VectorFunctionSpace, assign
    import sys

    # which forces you want to post-process: (subset of pnps.functionals and fkeys in simulation_script)
    fkeys = {"Javgtop","Javgctr","Javgbtm","Fp[0]","Fp[1]","Fp[2]","Fshear[0]",
            "Fshear[1]","Fshear[2]","Fbare[0]","Fbare[1]","Fbare[2]"}
    
    # --- get filenames and meshes
    
    # get git repo
    with open("%s/sim/gitrepo"%(DATADIR,), "r") as f:
        gitrepo = f.next().strip()
    
    # get data filename
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

    # get forces data
    force_iter = {}
    for k in fkeys:
        force_iter[k] = iter(data.data[k])

    # --- assign forces to Functions
    
    # initialize
    forces = {}
    D = FunctionSpace(mesh, "DG", 0)
    Ddofmap = D.dofmap()
    
    for k in fkeys:
        forces[k] = Function(D)
        
    # read from data to DG Functions
    for c in cells(mesh):
        i = Ddofmap.cell_dofs(c.index())
        for k in forces:
            forces[k].vector()[i] = force_iter[k].next()
            
    # interpolate to CG Functions
    V = FunctionSpace(mesh_o, "CG", 1)
    for k in forces:
        f = Function(V)
        f.interpolate(forces[k])
        forces[k] = f
        
    for i in range(3):
        f = Function(V)
        f.assign(sum(forces["%s[%d]"%(s,i)] for s in ["Fp", "Fshear", "Fbare"]))
        forces["F[%d]"%i] = f
        
    vforces = {}
    VV = VectorFunctionSpace(mesh_o, "CG", 1)
    
    for k in {"Fp", "Fshear", "Fbare", "F"}:
        forces[k] = [forces["%s[%d]"%(k,j)] for j in range(3)]
        vforces[k] = Function(VV)
        assign(vforces[k], [forces[k+"[0]"], forces[k+"[2]"]])
    
    #fmaps = {}
    #def make_force_map(k):
    #    return lambda x: map(lambda j : forces[k][j](x), [0,1,2])    
    #fmaps = {k:make_force_map(k) for k in {"Fp", "Fshear", "Fbare", "F"}}
    
    F = lambda x: map(lambda j : forces["F"][j](x), [0,1,2])
    
    if visualize:
        x = cells(mesh).next().midpoint()
        s = "F at %s: %s"
        print s % ([x[0],x[1]], F(x))
        print s % ([0.0,0.0], F([0.0, 0.0]))
      
        for k in vforces:
            plot(vforces[k], title=k)
 
        plot(forces["F[2]"], title="F[2]")
        interactive()
       
    return F

if __name__ == "__main__":
    F = post_process(visualize=True)
else:
    F = post_process()
