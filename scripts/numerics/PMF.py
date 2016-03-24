"plot 1D force/PMF profiles for 2D Howorka pore and save"

import os, numpy, dolfin, Howorka
from nanopores import kB, T, add_params, save_dict, saveplots
from matplotlib.pyplot import figure, plot, legend, show, title, xlabel, ylabel, savefig

add_params(
himp = .2,
hexp = .5,
Nimp = 1e5,
Nexp = 2e4,
z = (8., -8., 2)
)

# get force from explicit molecule
def F_explicit(*lspace):
    for z0 in numpy.linspace(*lspace):
        geo, phys = Howorka.setup2D(z0=z0, h=hexp)
        dolfin.plot(geo.boundaries, key="b", title="boundaries")
        pb, pnps = Howorka.solve2D(geo, phys, Nmax=Nexp, cheapest=True)
        yield pnps.zforces()
        
# get force from implicit molecule
def F_implicit(*lspace):
    geo, phys = Howorka.setup2D(z0=None, h=himp)
    pb, pnps = Howorka.solve2D(geo, phys, Nmax=Nimp, cheapest=True)
    for z0 in numpy.linspace(*lspace):
        yield pnps.zforces_implicit(z0)
    #pnps.visualize()
           
# compute PMF from force (antiderivative with trapezoid rule)
def PMF(F, a, b, N):
    dx = (b - a)/(N - 1.) * 1e-9 * 1e-12/(kB*T)
    forces = F(a, b, N)
    u = [0., 0., 0.]
    f0 = forces.next()
    yield tuple(u), f0
    for f1 in forces:
        for i in (0, 1, 2):
            u[i] -= 0.5*dx*(f0[i] + f1[i])
        f0 = f1
        yield tuple(u), f0
        
# plot force from implicit molecule
space = (9., -9., 101)
x = numpy.linspace(*space)
y = list(PMF(F_implicit, *space))
(u, uel, udrag), (f, fel, fdrag) = (
        tuple(tuple([yy[j][i] for yy in y] for i in (0, 1, 2)) for j in (0, 1)))

for i, ff in enumerate([u, uel, udrag, f, fel, fdrag]):
    figure(i)
    plot(x, ff, "-", label="implicit")

# plot force from explicit molecule and save figures
space = z
x = numpy.linspace(*space)
y = list(PMF(F_explicit, *space))

(u, uel, udrag), (f, fel, fdrag) = (
        tuple(tuple([yy[j][i] for yy in y] for i in (0, 1, 2)) for j in (0, 1)))

style = "s--"
label = "explicit"
xlab = "z-coordinate of molecule center [nm]"
ylabu = "PMF [kT]"
ylabf = "force [pN]"
folder = os.path.expanduser("~") + "/papers/pnps-numerics/figure_material/PMF/"
fname = folder + "%s.eps"
ifig = 0

def plotF(x, u, name, titl, ylab):
    global ifig
    figure(ifig)
    plot(x, u, style, label=label)
    title(titl)
    xlabel(xlab)
    ylabel(ylab)
    legend()
    savefig(fname % name, bbox_inches='tight')
    ifig += 1

plotF(x, u, "PMF", "PMF", ylabu)
plotF(x, uel, "PMFel", "PMF from electric force", ylabu)
plotF(x, udrag, "PMFdrag", "PMF from drag force", ylabu)
plotF(x, f, "F", "F", ylabf)
plotF(x, fel, "Fel", "Fel", ylabf)
plotF(x, fdrag, "Fdrag", "Fdrag", ylabf)

# save metadata
save_dict(PARAMS, folder, "meta")

# save plots
saveplots("HoworkaPMF", meta=PARAMS)

