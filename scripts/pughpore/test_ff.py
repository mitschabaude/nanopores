# (c) 2016 Gregor Mitscha-Baude
import numpy as np
import nanopores.models.pughpore as pugh

params = dict(h=2., Nmax=6e5)
H = pugh.pughpore.params["H"]
eps = 5.

ran = np.linspace(-H/2.+eps, H/2.-eps, 24)
X = [[0.,0.,t] for t in ran]

pugh.F_explicit(X, nproc=6, name="pughcenter", **params)
