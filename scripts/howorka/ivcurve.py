# iv curve or other stuff for Howorka where geometry can be reused
import numpy
import nanopores
import matplotlib.pyplot as plt
from nanopores.models import Howorka

nanopores.add_params(**Howorka.PARAMS)
nanopores.add_params(
    bV = numpy.linspace(-0.1,0.1,11),
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
    J = pnps.get_functionals(["Javgctr"])["Javgctr"]
    return dict(J=J)

#result, stamp = nanopores.iterate_in_parallel(run, iterkeys=[plot], **PARAMS)
plots = nanopores.parallel_output(run, showplot=False, **PARAMS)

ax = plots["J"]
line = ax.lines[0]
line.set_label("PNPS simulation")
ax.set_ylabel("current [pA]")
ax.set_xlabel("voltage [V]")

# plot experimental data from text file
csvfile = 'burns16iv.csv'
IV = numpy.genfromtxt(csvfile, delimiter=',')
mV = IV[:,0]
V = 1e-3*mV
I = IV[:,1]

ax.plot(V, I, "s--", label="Burns et al. '16")
ax.legend(loc="best")

fig = plt.gcf()
fig.set_size_inches(6, 6)
plt.show()
#nanopores.showplots()