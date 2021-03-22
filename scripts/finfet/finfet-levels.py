import dolfin, nanopores
from dolfin import *
from random import random
from nanopores.geometries.finfet import finfet, dopants
from nanopores import eperm, cm, qq, showplots, Estimator
from matplotlib.pyplot import loglog

Ndop = 12

def solve(h):
    # --- create mesh and geometrical/physical context information
    tic()
    t = dolfin.Timer("mesh")
    geo = finfet.create_geometry(lc=h)
    print("Number of elements:", geo.mesh.num_cells())
    print("Number of vertices:", geo.mesh.num_vertices())
    #finfet.plot()
    phys = nanopores.Physics("finfet", geo, dopants=dopants(Ndop), vD=None, vG=0., vS=None)
    #phys.add_dopants
    t.stop()

    # --- definition and solution of PDE
    t = dolfin.Timer("init")
    pde = nanopores.NonstandardPB(geo, phys)
    pde.tolnewton = 1e-5
    pde.newtondamp = 1.
    t.stop()
    t = dolfin.Timer("solve")
    pde.solve()
    t.stop()
    u = pde.solution
    return u, toc()
    
def hplot(name, h, err, xlab="", ylab="", rate=None, fig=True):
    from matplotlib.pyplot import figure, loglog, xlabel, ylabel, legend, show
    if fig is True:
        figure()
    loglog(h, err, 's-', label=name)
    if rate and h[0] != 0:
        alg = [err[0]/(h[0]**rate)*n**rate for n in h]
        loglog(h, alg, '--', label="h^(%s)" %rate)
    #xlabel("# Elements")
    xlabel(xlab)
    ylabel(ylab)
    legend(loc='upper right')

# --- solve for different h ---
solutions = []
times = []
hs = 2., 1.5, 1.2, 1., .8, .6, .5, .25

print("---- SOLVING ----")
for h in hs:
    print("-- h = %s --" %h)
    u, time = solve(h)
    print("\nTIME:", time, "s\n")
    solutions.append(u)
    times.append(time)
    dolfin.plot(u, title="potential, h=%s" %h)

print("---- EXTRAPOLATING SOLUTION AND COMPUTING NORMS ----")
# with extrapolation
"""
tic()
ulast = solutions[-1]
mesh = ulast.function_space().mesh()
V = FunctionSpace(mesh, "CG", 2)
u = Function(V)
u.extrapolate(ulast)
deg = 1
print "Time for extrapolation:", toc(), "s"
"""
# without extrapolation
ulast = solutions.pop()
tlast = times.pop()
hs = hs[:-1]
deg = 0
print("Time for last solution:", tlast, "s")

def err(uh):
    return errornorm(u, uh, "H1", degree_rise=deg)
tic()
errors = [err(uh) for uh in solutions]
print("Time for norms:", toc(), "s")

hplot("error", hs, errors, xlab="h [nm]", ylab="H1 error", rate=1.)
hplot("work", hs, times, xlab="h [nm]", ylab="work [s]", rate=-3.)

print(times, errors)

showplots()
list_timings()
interactive() 

