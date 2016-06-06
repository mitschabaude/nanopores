"""
TODO:
    -) understand unit in which current is measured
    -) plot current for different radii and Qs
"""
import nanopores, dolfin, Howorka
from nanopores.physics.simplepnps import SimpleNernstPlanckProblem
import matplotlib.pyplot as plt
import force_profiles
#from force_profiles import geo, loadall, function_from_lambda

nanopores.add_params(
rMol = 0.5,
)

name = {0.5: "r05", 0.2: "r02"}[rMol]

class DiffusionProblem1D(SimpleNernstPlanckProblem):
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
        
        c = dolfin.TrialFunction(V)
        d = dolfin.TestFunction(V)
        
        FF = dolfin.as_vector([F])
        J = -D*grad(c) + D/kT*FF*c
        a = dolfin.inner(J, grad(d))*dx
        L = dolfin.Constant(0.)*d*dx
        return a, L
            
    @staticmethod
    def bcs(V, geo, c0):
        return geo.pwBC(V, "c0", value={"top" : c0, "bottom": c0})
        
def current_per_ns(geo, phys, c, F):
    dx = geo.dx()
    grad = phys.grad
    kT = dolfin.Constant(phys.kT)
    D = dolfin.Constant(Dtarget(phys.rTarget))
    FF = dolfin.as_vector([F])
    print "v = %s" % (Dtarget(phys.rTarget)*F(0.)/phys.kT,)
    
    j = (-D*grad(c) + D/kT*FF*c)
    #dolfin.plot(j)
    #dolfin.interactive()
    L = 20.
    J = dolfin.assemble(j[0]/dolfin.Constant(L) * geo.dx())
    return J
        
def Dtarget(r):
    return nanopores.kT/(6*dolfin.pi*nanopores.eta*r)
    
def J_FEM(F, c0):
    geo = force_profiles.geo
    phys = nanopores.Physics(geo=geo, rTarget=rMol*1e-9, lscale=1e9)
    pde = nanopores.solve_pde(DiffusionProblem1D, geo=geo, phys=phys, F=F, c0=c0, verbose=False)
    c = pde.solution
    return c, current_per_ns(geo, phys, c, F)

c0 = 1.
for F, Fi, Fi_better, alpha, beta, Q in force_profiles.Forces(name):
    F = force_profiles.function_from_lambda(lambda z: 1e-12*F(z))
    u, J = J_FEM(F, c0)
    
    
    F = force_profiles.function_from_lambda(lambda z: 1e-12*Fi(z))
    ui, Ji = J_FEM(F, c0)  
    
    F = force_profiles.function_from_lambda(lambda z: 1e-12*Fi_better(z))
    uib, Jib = J_FEM(F, c0)
    
    print "Q %s, J %s, Ji %s, Jib %s" % (Q, J, Ji, Jib)
    
    #force_profiles.plot_function(F, label="Q="+str(Q))
    force_profiles.plot_function(u, label="Q="+str(Q))

plt.legend()
plt.show()




