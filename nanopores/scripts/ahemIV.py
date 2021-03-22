import math, nanopores, dolfin

# @Benjamin, Gregor TODO:
# -) check permittivity and surface charge of ahem
# -) what biased voltage to use?

# some default values for parameters

### geo params [nm]
geo_name = "aHem"
domscale = 1.
l4 = 15.
l3 = 15.
R = 20.
r0 = 5.
z0 = 10.
exit_i = 1
badexit = {"upperbulkb"}
goodexit = {"exit"}

### phys params
phys_name = "pore_molecule"
bV = .5 # [V]
ahemqs = 0.0 # [C/m**2]
rTarget = 0.5e-9 # [m] for diffusion coeff.
bulkcon = 1000

### num params
clscale = 10.
refinement = False
maxcells = 50e3
newtondamp = 1.0
reuse_mesh = False
tolnewton = 1e-1
skip_stokes = True
iterative = True

def _update(dic, dic2): # conservative update
    dic.update({key:dic2[key] for key in dic2 if not key in dic})

def _globals(): # globals except hidden ones ("_"), modules and functions
    from types import ModuleType, FunctionType
    return {key : var for key, var in list(globals().items()) if not
        (key.startswith("_") or isinstance(var, ModuleType) or isinstance(var, FunctionType))}

def calculate(**params):
    # this time we do it the simple way: just pass every parameter
    # have to be careful though
    globals().update(params)
    params.update(_globals())

    # use some of the parameters
    params["x0"] = [r0, 0., z0]
    params["l3"] = l3*domscale
    params["R"] = R*domscale
    # TODO does this something?
    nanopores.IllposedNonlinearSolver.newtondamp = newtondamp
    nanopores.PNPS.tolnewton = tolnewton

    t = dolfin.Timer("meshing")
    geo = nanopores.geo_from_xml_threadsafe(geo_name, **params)
    print("Mesh generation time:",t.stop())

    #dolfin.plot(geo.submesh("solid"), interactive=True)
    phys = nanopores.Physics(phys_name, geo, **params)

    t = dolfin.Timer("PNPS")
    pnps = nanopores.PNPS(geo, phys)
    if skip_stokes:
        pnps.solvers.pop("Stokes")
    pnps.alwaysstokes = True
    pnps.solve()
    print("Time to calculate F:",t.stop())
    #pnps.visualize("fluid")

    (v, cp, cm, u, p) = pnps.solutions(deepcopy=True)
    # F = phys.Feff(v, u)
    # def avg(u, dx):
    #     return dolfin.assemble(u*dx)/dolfin.assemble(dolfin.Constant(1.)*dx)

    Jcomp = ["Jzdiff", "Jzdrift", "Jzstokes"]
    lPore = geo.params["ltop"]+geo.params["lctr"]+geo.params["lbtm"]
    Jzdiff = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.D*phys.rDPore*phys.grad(-cp+cm)[2] /lPore * geo.dx("pore")
    Jzdrift = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.mu*phys.rDPore*(-cp-cm)*phys.grad(v)[2]/lPore * geo.dx("pore")
    Jzstokes = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.stokesdampPore*(cp-cm)*u[2]/lPore * geo.dx("pore")

    Jcomponents = [j+p for j in Jcomp for p in ["top","btm"]]
    Jzdifftop = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.D*phys.rDPore*phys.grad(-cp+cm)[2] /geo.params["ltop"] * geo.dx("poretop")
    Jzdrifttop = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.mu*phys.rDPore*(-cp-cm)*phys.grad(v)[2]/geo.params["ltop"] * geo.dx("poretop")
    Jzstokestop = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.stokesdampPore*(cp-cm)*u[2]/geo.params["ltop"] * geo.dx("poretop")
    Jzdiffbtm = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.D*phys.rDPore*phys.grad(-cp+cm)[2] /geo.params["lbtm"] * geo.dx("porebottom")
    Jzdriftbtm = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.mu*phys.rDPore*(-cp-cm)*phys.grad(v)[2]/geo.params["lbtm"] * geo.dx("porebottom")
    Jzstokesbtm = dolfin.Constant((1.0/phys.lscale)**2) * phys.cFarad*phys.stokesdampPore*(cp-cm)*u[2]/geo.params["lbtm"] * geo.dx("porebottom")

    result = pnps.get_functionals()

    for j in Jcomp+Jcomponents:
        result.update({j: 1e12*dolfin.assemble(locals()[j])})

    return result
