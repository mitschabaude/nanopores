# (c) 2016 Gregor Mitscha-Baude
import nanopores.models.pughpore as pugh

params = dict(h=2., Nmax=2e2)
X = [[0.,0.,float(t)] for t in [23.]]

print pugh.F_explicit(X, nproc=1, name="testpugh", **params)
