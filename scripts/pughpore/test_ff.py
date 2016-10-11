import nanopores.models.pughpore as pugh
from nanopores.tools.solvers import calculate_forcefield

calculate = pugh.F_explicit
params = dict(h=10., Nmax=1e2)
default = dict(pugh.default.geop)
default.update(pugh.default.physp)
default.update(pugh.default.solverp)

name = "testpugh"
X = [[0.,0.,t] for t in [1.,2.]]
calculate_forcefield(name, X, calculate, params, default)
