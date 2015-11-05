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
r0 = 0.
z0 = 10.
exit_i = None

### phys params
phys_name = "pore_molecule"
bV = .5 # [V]
ahemqs = 0.0 # [C/m**2]
rTarget = 0.5e-9 # [m] for diffusion coeff.
bulkcon = 3e2

### num params
clscale = 10.
refinement = False
maxcells = 50e3
newtondamp = 1.0
reuse_mesh = False
tolnewton = 1e-1

def _update(dic, dic2): # conservative update
    dic.update({key:dic2[key] for key in dic2 if not key in dic})
    
def _globals(): # globals except hidden ones ("_"), modules and functions
    from types import ModuleType, FunctionType
    return {key : var for key, var in globals().items() if not
        (key.startswith("_") or isinstance(var, ModuleType) or isinstance(var, FunctionType))}

def calculate(**params):
    # this time we do it the simple way: just pass every parameter
    # have to be careful though
    _update(params, _globals())
    globals().update(params)
    
    # use some of the parameters
    params["x0"] = [r0, 0., z0]
    params["l3"] = l3*domscale
    params["R"] = R*domscale
    nanopores.IllposedNonlinearSolver.newtondamp = newtondamp
    nanopores.PNPS.tolnewton = tolnewton
    
    t = dolfin.Timer("meshing")
    geo = nanopores.geo_from_xml_threadsafe(geo_name, **params)
    print "Mesh generation time:",t.stop()
    
    #dolfin.plot(geo.submesh("solid"), interactive=True)
    phys = nanopores.Physics(phys_name, geo, **params)

    #print "Physics:"
    #for item in phys.__dict__.items():
    #    print "%s = %s" %item
    
    pnps = nanopores.PNPS(geo, phys)
    pnps.solve()
    pnps.visualize("fluid")
    
    (v, cp, cm, u, p) = pnps.solutions(deepcopy=True)
    F = phys.Feff(v, u)

    # some points to evaluate exit time
    x0 = geo.params["x0"]
    rx0 = math.sqrt(sum(x**2 for x in x0))
    rnear = rx0 - geo.params["rMolecule"]
    rfar = rx0 + geo.params["rMolecule"]
    xnear = map(lambda x: rnear/rx0*x, x0)
    xfar = map(lambda x: rfar/rx0*x, x0)
    xtop = [0.,0.,0.]
    xbtm = [0.,0.,geo.params["zbtm"]]

    def avg(u, dx):
        return dolfin.assemble(u*dx)/dolfin.assemble(dolfin.Constant(1.)*dx)

    etp = nanopores.LinearPDE(geo, nanopores.ExitTimeProblem, phys, F=F)
    etp.solve()
    tau = etp.solution
    
    result = dict(
        tmolmin = tau(xnear),
        tmolavg = avg(tau, geo.dS("moleculeb")),
        tmolmax = tau(xfar),
        tporemax = tau(xtop),
        tporemin = tau(xbtm),
        tporeavg = avg(tau, geo.dx("pore")),
    )
    return result


