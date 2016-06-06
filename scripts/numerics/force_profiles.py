"calculate and save 1D force profiles for 2D Howorka pore for different molecule charges"

import os, numpy, dolfin
import nanopores, Howorka
import matplotlib.pyplot as plt
#from nanopores import kB, T, add_params, save_dict, saveplots, showplots
#from matplotlib.pyplot import figure, plot, legend, show, title, xlabel, ylabel, savefig

nanopores.add_params(
himp = .2,
hexp = .5,
Nimp = 5e4,
Nexp = 2e4,
Nz = 50,
name = "",
save = False,
)

Qmols = [-2., -1., 1e-4, 1., 2.]
a, b = -10., 10.
folder = os.path.expanduser("~") + "/papers/pnps-numerics/data/forces/"

space = (a, b, Nz)

domain = nanopores.Interval(a, b)
domain.addsubdomains(fluid = domain)
domain.addboundaries(top = domain.boundary("right"),
                     bottom = domain.boundary("left"))
                     
geo = domain.create_geometry(lc=(b-a)/(Nz-1))
mesh = geo.mesh
V = dolfin.FunctionSpace(mesh, "CG", 1)
v2d = dolfin.vertex_to_dof_map(V)
coord = geo.mesh.coordinates()[:, 0] # numpy array
linspace = numpy.linspace(*space)

def function_from_values(values):
    u = dolfin.Function(V)
    u.vector()[v2d] = numpy.array(values)
    return u
    
def function_from_lambda(f):
    values = [f(z) for z in coord]
    return function_from_values(values)
    
def plot_function(u, *args, **kwargs):
    plt.plot(linspace, [u(z) for z in linspace], *args, **kwargs)
    
def plot_point(F):
    plot_function(F, "-", label="point-sized")
def plot_finite(F):
    plot_function(F, "s--", label="finite-sized")    
def post_plot():
    plt.xlabel("z-coordinate of molecule center [nm]")
    plt.ylabel("force [pN]")
    plt.legend(loc="best")

def plot_profile(Fimp, Fexp):
    plt.figure()
    plot_point(Fimp)
    plot_finite(Fexp)
    post_plot()
    
# TEST:
#plot_function(function_from_lambda(numpy.sin))
#nanopores.showplots()
#exit()

# get force from explicit molecule
def F_explicit(Qmol):
    values = []
    for z0 in coord:
        geo, phys = Howorka.setup2D(z0=z0, h=hexp, Qmol=Qmol)
        dolfin.plot(geo.boundaries, key="b", title="boundaries")
        pb, pnps = Howorka.solve2D(geo, phys, Nmax=Nexp, cheapest=True)
        dolfin.plot(geo.boundaries, key="b", title="boundaries")
        values.append(pnps.zforces())
    F, Fel, Fdrag = tuple(function_from_values(v) for v in zip(*values))
    return F, Fel, Fdrag
        
# get force from implicit molecule
def F_implicit(Qmol):
    geo, phys = Howorka.setup2D(z0=None, h=himp, Qmol=Qmol)
    pb, pnps = Howorka.solve2D(geo, phys, Nmax=Nimp, cheapest=True)
    values = [pnps.zforces_implicit(z) for z in coord]
    F, Fel, Fdrag = tuple(function_from_values(v) for v in zip(*values))
    #pnps.visualize()
    return F, Fel, Fdrag
    
def saveforces(name, F, Fel, Fdrag):
    #dolfin.File(name + "_mesh.xml") << geo.mesh
    dolfin.File(folder+name + "_F.xml") << F
    dolfin.File(folder+name + "_Fel.xml") << Fel
    dolfin.File(folder+name + "_Fdrag.xml") << Fdrag
    
def loadforces(name):
    #mesh = Mesh(name + "_mesh.xml")
    #V = FunctionSpace(mesh, "CG", 1)
    F = dolfin.Function(V, folder+name + "_F.xml")
    Fel = dolfin.Function(V, folder+name + "_Fel.xml")
    Fdrag = dolfin.Function(V, folder+name + "_Fdrag.xml")
    return F, Fel, Fdrag
    
def saveall(name, Qmol):
    name = name + "_Q%.2f" % Qmol
    Fi, Feli, Fdragi = F_implicit(Qmol=Qmol)
    saveforces(name + "_imp", Fi, Feli, Fdragi)
    F, Fel, Fdrag = F_explicit(Qmol=Qmol)
    saveforces(name + "_exp", F, Fel, Fdrag)
    
def loadall(name, Qmol):
    name = name + "_Q%.2f" % Qmol
    imp = loadforces(name + "_imp")
    exp = loadforces(name + "_exp")
    return imp, exp

if save:
    for Q in Qmols:
        saveall(name, Q)
    exit()
    
def construct_alpha(a0):
    lpore = 4.5 # TODO
    a = 4.
    b = 6.
    def alpha(z):
        if abs(z) > b:
            return 1.
        elif abs(z) < a:
            return a0
        else: # abs(z) in [a,b]
            return 1. + (b-abs(z))/(b-a)*(a0 - 1)
            
    return function_from_lambda(alpha)
    
def Forces(name):
    for Q in Qmols:
        (Fi, Feli, Fdragi), (F, Fel, Fdrag) = loadall(name, Q)
        alpha0 = Fdrag(0.0)/Fdragi(0.0)
        alpha = construct_alpha(alpha0)
        beta0 = Fel(0.0)/Feli(0.0)
        beta = construct_alpha(beta0)
        
        Fdragi_better = function_from_lambda(lambda z : Fdragi(z)*alpha(z))
        Feli_better = function_from_lambda(lambda z : Feli(z)*beta(z))
        Fi_better = function_from_lambda(lambda z : Feli_better(z) + Fdragi_better(z))
        
        yield F, Fi, Fi_better, alpha, beta, Q
        
if __name__ == "__main__":
    for F, Fi, Fi_better, alpha, beta, Q in Forces(name):        
        print "Q %s, alpha %s" % (Q, alpha(0.))
        plt.figure(0)
        plot_function(alpha, label="Q = %.0f"%Q)
        
        plt.figure(1)
        plot_function(beta, label="Q = %.0f"%Q)
        
        plt.figure()
        plot_finite(F)
        plot_point(Fi)
        plot_function(Fi_better, "-", label="point-sized, corrected")
        post_plot()
        
    plt.figure(0)
    plt.legend()
    plt.figure(1)
    plt.legend()
    nanopores.showplots()

