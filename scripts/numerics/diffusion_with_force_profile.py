"""
TODO:
    -) understand no boundary condition
    -) validate understanding with analytical solution
"""
import nanopores, dolfin, os
from nanopores.physics.simplepnps import SimpleNernstPlanckProblem
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import force_profiles
from collections import defaultdict

nanopores.add_params(
savefig = False
)

class DiffusionProblem1D(SimpleNernstPlanckProblem):
    method = SimpleNernstPlanckProblem.method
    method["iterative"] = False
    
    @staticmethod
    def initial_u(V, c0):
        u = dolfin.Function(V)
        u.interpolate(dolfin.Constant(c0))
        return u

    @staticmethod
    def forms(V, geo, phys, F):
        dx = geo.dx()
        grad = phys.grad
        kT = dolfin.Constant(phys.kT)
        D = dolfin.Constant(Dtarget(phys.rTarget))
        lscale = dolfin.Constant(phys.lscale)
        n = dolfin.FacetNormal(geo.mesh)
        
        c = dolfin.TrialFunction(V)
        d = dolfin.TestFunction(V)
        
        FF = dolfin.as_vector([F])
        J = -D*grad(c) + D/kT*FF*c
        a = dolfin.inner(J, grad(d))*dx
        L = dolfin.Constant(0.)*d*dx
        
        aNoBC = -lscale*dolfin.inner(J, n*d)*geo.ds("bottom")
        a += aNoBC
        return a, L
            
    @staticmethod
    def bcs(V, geo, c0):
        bc = dict(
            top = c0,
            #bottom = c0,
        )
        return geo.pwBC(V, "c0", value=bc)
        
def current(geo, phys, c, F):
    dx = geo.dx()
    grad = phys.grad
    lscale = phys.lscale
    mol = phys.mol
    kT = dolfin.Constant(phys.kT)
    D = dolfin.Constant(Dtarget(phys.rTarget))
    FF = dolfin.as_vector([F])
    print "v = %s" % (Dtarget(phys.rTarget)*F(0.)/phys.kT,)
    
    j = -D*grad(c) + D/kT*FF*c
    #dolfin.plot(j)
    #dolfin.interactive()
    L = 20.
    r0 = 1./lscale
    Across = r0**2 * dolfin.pi
    # current in N/s
    J = mol * Across * dolfin.assemble(j[0]/dolfin.Constant(L) * dx)
    # current in N/ms
    J = J * 1e-3
    return J
        
def Dtarget(r):
    return nanopores.kT/(6*dolfin.pi*nanopores.eta*r)
    
def J_FEM(F, c0):
    geo = force_profiles.geo
    phys = nanopores.Physics(geo=geo, rTarget=rMol*1e-9, lscale=1e9)
    pde = nanopores.solve_pde(DiffusionProblem1D, geo=geo, phys=phys,
                              F=F, c0=c0, verbose=False)
    c = pde.solution"Simulation"+
    return c, current(geo, phys, c, F)

def gather_currents(name, c0):
    currents = defaultdict(list)
    qmols = []
    
    for results in force_profiles.Forces(name):
        qmols.append(results["Q"])    
        
        for key in "F", "Fi", "Fi2":
            f = results[key]
            f = force_profiles.function_from_lambda(lambda z: 1e-12*f(z))
            u, J = J_FEM(f, c0)
            currents[key].append(J)
            #force_profiles.plot_function(f, label="Q="+str(Q))
            #if key=="F":
            #    force_profiles.plot_function(u, label="Q="+str(results["Q"]))
        
        print "Q %s, J %s, Ji %s, Jib %s" % (
        qmols[-1], currents["F"][-1], currents["Fi"][-1], currents["Fi2"][-1])
        
    return qmols, currents
    
c0 = 1.6605 # [mol/m**3] = 1 molecule per (10nm)**3
#names = {0.25: "r025", 0.5: "r05", 0.2: "r02", 0.4: "r04", 0.75: "r075"}
items = (0.25, "r025"), (0.5, "r05"), (0.75, "r075")
    
figures = os.path.expanduser("~") + "/papers/pnps-numerics/figures/"    

    
for rMol, name in items:
    #plt.figure()
    qmols, currents = gather_currents(name, c0)
    #plt.legend()
    
    fig, ax = plt.subplots()    
    ax.plot(qmols, currents["F"], "s-g", label="finite")
    ax.plot(qmols, currents["Fi"], "s-b", label="point")
    ax.plot(qmols, currents["Fi2"], "s-r", label="point, corrected")

    #ax.set_aspect('equal')
    tick_spacing = 1.
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
    #ax.set_ylim([-0.4, 0.]) #xaxis.set_ticks(np.arange(start, end, 0.712123))
    #ax.grid(True, which='both')
    #ax.axhline(y=0, color='#cccccc')
      
    plt.title("r = %s" %rMol)
    plt.xlabel("Molecule charge [q]")
    plt.ylabel("Molecule current [1/ms]")
    plt.legend()
    if savefig:
        fig = plt.gcf()
        fig.set_size_inches(5, 4)
        plt.savefig(figures + "molcurrent_r%.2f.eps" % rMol, bbox_inches='tight')
    
#plt.show()




