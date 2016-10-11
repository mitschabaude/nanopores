# (c) 2016 Gregor Mitscha-Baude
"PNPS solvers and visualization for pugh pore"

import dolfin
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
import nanopores.physics.simplepnps as simplepnps
import nanopores.tools.box as box
from nanopores.tools.solvers import Setup as SetupBase

default = nano.Params(
geop = nano.Params(
    R = 30.,
    H = 80.,
    x0 = [0.,0.,15.],
),
physp = nano.Params(
    Qmol = -1.,
    bulkcon = 300.,
    dnaqsdamp = .5,
    bV = -0.1,
),
solverp = nano.Params(
    h = 2.,
    frac = 0.2,
    Nmax = 5e5,  
    imax = 30,
    tol = 1e-2,
    cheapest = False,
))
# change global settings of mesh generation
box.set_tol(None)

class Setup(SetupBase):
    default = default
                    
    def init_geo(self):
        geo = pughpore.get_geo(self.solverp.h, **self.geop)
        molec = nano.curved.Sphere(geo.params["rMolecule"], geo.params["x0"])
        geo.curved = dict(moleculeb = molec.snap)
        self.geo = geo
        
    def init_phys(self):
        self.phys = nano.Physics("pore_mol", self.geo, **self.physp)
        set_sideBCs(self.phys, self.geop, self.physp)

def solve(setup, visualize=False):
    geo, phys, solverp = setup.geo, setup.phys, setup.solverp
    if visualize:
        R, H = geo.params["R"], geo.params["H"]
        mesh2D = nano.RectangleMesh([-R,-H/2.], [R, H/2.], int(4*R), int(2*H))
    else:
        mesh2D = None
    if geo.mesh.num_cells() < solverp.Nmax:
        pb = prerefine(setup, visualize, mesh2D)
    else:
        pb = None
        
    pnps = simplepnps.PNPSFixedPointbV(geo, phys, ipicard=solverp.imax,
               verbose=True, #taylorhood=True,
               stokesiter=True, tolnewton=solverp.tol, iterative=True)          
    
    print "Number of cells:", geo.mesh.num_cells()
    print "DOFs:", pnps.dofs()
    dolfin.tic()
    for i in pnps.fixedpoint(ipnp=5):
        if visualize:
            v, cp, cm, u, p = pnps.solutions()
            nano.plot_cross(v, mesh2D, title="potential", key="u")
            #nano.plot_cross(cm, mesh2D, title="negative ions", key="cm")
    print "CPU time (solve): %.3g s" %(dolfin.toc(),)
    #if visualize:
    #    nano.plot_cross_vector(u, mesh2D, title="u")
    return pb, pnps

def get_forces(setup, pnps):
    forces = pnps.evaluate(setup.phys.CurrentPNPS)
    forces.update(pnps.evaluate(setup.phys.ForcesPNPS))
    return forces
      
def prerefine(setup, visualize=False, mesh2D=None):
    geo, phys, p = setup.geo, setup.phys, setup.solverp
    dolfin.tic()
    goal = phys.CurrentPB
    pb = simplepnps.SimpleLinearPBGO(geo, phys, goal=goal, cheapest=p.cheapest)
    
    for i in pb.adaptive_loop(p.Nmax, p.frac):
        if visualize:
            dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
            #nano.plot_cross(pb.solution, mesh2D,title="pb potential", key="pb")
    print "CPU time (PB): %.3g s" %(dolfin.toc(),)
    return pb    

def solve1D(geop, physp):
    geo = pughpore.get_geo1D(lc=.001, **geop)
    phys = nano.Physics("pore", geo, **physp)
    pnp = nano.solve_pde(simplepnps.SimplePNPProblem, geo, phys)
    return geo, pnp
    
def visualize1D(geo, pnp):
    v, cp, cm = pnp.solutions()
    h = geo.params["H"]
    nano.plot1D({"potential": v}, (-h/2, h/2, 1001),
                "x", dim=1, axlabels=("z [nm]", "potential [V]"))
    nano.plot1D({"c+": cp, "c-":cm},  (-h/2, h/2, 1001),
                "x", dim=1, axlabels=("z [nm]", "concentrations [mol/m^3]"))
                
class u1D(dolfin.Expression):
    def __init__(self, u, damping=1.):
        self.u = u
        self.damping = damping
        dolfin.Expression.__init__(self)
        
    def damp(self, scalar):
        self.damping *= scalar

    def eval(self, value, x):
        value[0] = self.damping*self.u(x[2])
        
def set_sideBCs(phys, geop, physp):
    geo, pnp = solve1D(geop, physp)
    v, cp, cm = pnp.solutions()
    phys.v0["sideb"] = u1D(v)
    phys.cp0["sideb"] = u1D(cp)
    phys.cm0["sideb"] = u1D(cm)
    
def join_dicts(list):
    # [{"F":1.0}, {"F":2.0}, ...] --> {"F":[1.0, 2.0, ...]}
    return {key:[dic[key] for dic in list] for key in list[0]}
    
# evaluate finite-size model for a single x0 position
def F_explicit(x0, **params):
    setup = Setup(x0=x0, **params)
    pb, pnps = solve(setup, True)
    return get_forces(setup, pnps)
    
     
# TODO
## evaluate point-size model for a number of z positions
#def F_implicit3D(Z, **params):
#    geo, phys = setup2D(z0=None, **params)
#    pb, pnps = solve2D(geo, phys, **params)
#    values = [pnps.zforces_implicit(z0) for z0 in Z]
#    F, Fel, Fdrag = tuple(zip(*values))
#    return F, Fel, Fdrag
#
## get discrete force fields from point-size model
#def F_field_implicit3D(**params):
#    params["z0"] = None
#    geo, phys = setup2D(**params)
#    pb, pnps = solve2D(geo, phys, **params)
#    (v, cp, cm, u, p) = pnps.solutions()
#    F, Fel, Fdrag = phys.Forces(v, u)
#    return F, Fel, Fdrag
    
if __name__ == "__main__":
    setup = Setup(h=10., Nmax=1e5, x0=None)
    _, pnps = solve(setup, True)
    print get_forces(setup, pnps)
    dolfin.interactive()
    