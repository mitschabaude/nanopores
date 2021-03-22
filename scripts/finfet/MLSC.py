import numpy, dolfin, nanopores
from dolfin import *
from nanopores.geometries.finfet import finfet
from nanopores import showplots, saveplots, add_params
from collocation import dopants

add_params(
Ndop = 4,
h = 1.,
maxorder = 2,
)

def solve(geo, dops):
    tic()
    t = dolfin.Timer("phys")
    phys = nanopores.Physics("finfet", geo,
        dopants=dops,
        vD=None, vG=None, vS=None)
    phys.add_dopants
    t.stop()
    dolfin.plot(geo.submesh("source"), key="dop", title="dopants in source")
    
    t = dolfin.Timer("init")
    pde = nanopores.NonstandardPB(geo, phys)
    pde.tolnewton = 1e-8
    pde.newtondamp = 1.
    t.stop()
    t = dolfin.Timer("solve")
    pde.solve()
    t.stop()
    u = pde.solution
    print("Loop time: %s s" %(toc(),))
    return u
    
def interpolate(h, order):
    print("\nMeshing.")
    geo = finfet.create_geometry(lc=h)
    print("Number of elements:", geo.mesh.num_cells())
    print("Number of vertices:", geo.mesh.num_vertices())
    dopants_, weights = dopants(Ndop, order)
    Nsamples = len(dopants_)
    V = FunctionSpace(geo.mesh, "CG", 1)
    Iluh = Function(V)
    for i, (dops, weight) in enumerate(zip(dopants_, weights)):
        print("\nSample %d of %d:\nDopants:\n" %(i+1, Nsamples), dops, "\n")
        geo = finfet.recreate_geometry()
        u = solve(geo, dops)
        Iluh.vector()[:] += weight * u.vector()[:]
    return Iluh, Nsamples
    
def hplot(name, h, err, xlab="", ylab="", hstr="h", rate=None, fit=True, fig=True):
    from matplotlib.pyplot import figure, loglog, xlabel, ylabel, legend, show
    if fig is True:
        figure()
    loglog(h, err, 's-', label=name)
    if fit:
        p = numpy.polyfit(numpy.log(h), numpy.log(err), 1)
        C = numpy.exp(p[1])
        alpha = p[0]
        alg = [C*n**alpha for n in h]
        loglog(h, alg, '--', label="%.3f*%s^(%.3f)" %(C, hstr, alpha))
    if rate and h[0] != 0:
        alg = [err[0]/(h[0]**rate)*n**rate for n in h]
        loglog(h, alg, '--', label="%s^(%.3f)" %(hstr, rate))
    xlabel(xlab)
    ylabel(ylab)
    legend(loc='best')

# --- create interpolation for different h and orders ---
orders = list(range(1, maxorder+1))
samples = []
interpolands = []

for order in orders:
    print("\n-- order = %s --" %order)
    u, N = interpolate(h=h, order=order)
    samples.append(N)
    interpolands.append(u)

dolfin.plot(u, key="u", title="interpolated potential, order %s" %order)
ulast = interpolands.pop()
Nlast = samples.pop()

def err(uh):
    return errornorm(ulast, uh, "H1", degree_rise=0)/norm(ulast, "H1")

errors = [err(uh) for uh in interpolands]

label = "|Iluh - I%duh|/|I%duh|" %(maxorder, maxorder)
hplot(label, samples, errors, xlab="# samples M", ylab="relative interpolation error", hstr="M", fit=True)

list_timings()
#saveplots("mlsc", meta=PARAMS)
showplots()
interactive()

