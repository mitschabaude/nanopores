# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from nanopores import user_param

N = user_param(N=10)

ran = np.linspace(-30, 30, N)
X = [[0.,0.,t] for t in ran]
for x in X:
    print("[%.2f,%.2f,%.2f]" % tuple(x), end=' ')