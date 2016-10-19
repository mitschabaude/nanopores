# (c) 2016 Gregor Mitscha-Baude
import numpy as np
import nanopores.models.pughpore as pugh

params = dict(h=2., Nmax=6e5)
#params = dict(h=2., Nmax=6e5)
H = pugh.pughpore.params["H"]
eps = 5.

ran = np.linspace(-H/2.+eps, H/2.-eps, 48)
X = [[0.,0.,t] for t in ran]

hpore = pugh.pughpore.params["hpore"]
ran2 = np.linspace(-hpore/2.-1., hpore/2.+4, 12)
X2 = [[0.,0.,t] for t in ran2]

pugh.F_explicit(X+X2, nproc=6, name="pughcenter", **params)
