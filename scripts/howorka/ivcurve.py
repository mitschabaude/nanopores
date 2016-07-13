# iv curve or other stuff for Howorka where geometry can be reused
import numpy
import nanopores
import matplotlib.pyplot as plt
from nanopores.models import Howorka
from matplotlib2tikz import save as tikz_save

nanopores.add_params(**Howorka.PARAMS)
nanopores.add_params(
    bV = numpy.linspace(-0.1,0.1,11),
    dnaqsdamp = [0.2, 0.35, 0.5],
    bulkcon = 300,
    plot = "bV",
    nproc = 4,
)
print PARAMS

geo_params = dict(z0=None, rMolecule=rMolecule, Rx=Rx, Ry=Ry)
geo, phys = Howorka.setup2D(**geo_params)
mesh = geo.mesh

def run(**phys_params):
    params = geo_params.copy()
    params.update(phys_params)
    geo, phys = Howorka.setup2D(mesh=mesh, **params)
    pb, pnps = Howorka.solve2D(geo, phys, **params)    
    return dict(J = pnps.get_functional("Javgctr"))

#result, stamp = nanopores.iterate_in_parallel(run, iterkeys=[plot], **PARAMS)
plots = nanopores.parallel_output(run, showplot=False, **PARAMS)

# modify plot output
ax = plots["J"]
for line in ax.lines:
    label = line.get_label()
    #label = r"Sim." + (r"" if label.startswith("_") else r", %s" % label)
    label = label.replace("dnaqsdamp=", "$\\rho=-")
    label += "$" #r"q/\rm{nm}^2$"
    line.set_label(label)
ax.set_ylabel("current [pA]")
ax.set_xlabel("voltage [V]")

# plot experimental data from text file
csvfile = 'burns16iv.csv'
IV = numpy.genfromtxt(csvfile, delimiter=',')
mV = IV[:,0]
V = 1e-3*mV
I = IV[:,1]

ax.plot(V, I, "s", label="Burns et al.")
ax.set_xlim([bV[0]-0.01, bV[-1]+0.01])
ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#ax.legend(loc="best")

fig = plt.gcf()
fig.set_size_inches(2.6, 2.2)

# save to paper directory
from folders import FIGDIR
#tikz_save(FIGDIR+"iv2.tex", figureheight='2.5in', figurewidth='2.6in')
#plt.savefig(FIGDIR + "iv.eps", bbox_inches='tight')
plt.show()