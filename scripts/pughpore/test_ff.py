# (c) 2016 Gregor Mitscha-Baude
import nanopores.models.pughpore as pugh
from nanopores.tools.solvers import calculate_forcefield

calculate = pugh.F_explicit
params = dict(h=10., Nmax=2e2, dnaqsdamp=0.05)
default = dict(pugh.default.geop)
default.update(pugh.default.physp)
default.update(pugh.default.solverp)

name = "testpugh"
X = [[0.,0.,float(t)] for t in range(12)]
calculate_forcefield(name, X, calculate, params, default, nproc=6)
