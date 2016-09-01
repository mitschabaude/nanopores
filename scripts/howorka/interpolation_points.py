"create 2D point set for Howorka model where force shall be evaluated."
import numpy as np
from itertools import product
import math
gauss = np.polynomial.legendre.leggauss

# parameters TODO
h = 0.2
hout = 0.3

r0 = 1.
r = 0.5
l0 = 4.5
Ry = 7.
Rx = 3.

# effective pore radius for molecule
R = r0 - r

# get scaled 1D gauss quadrature nodes
# for x grid inside pore
# this has h = h/2 to resolve the thin region
k = int(math.ceil(R/h*2))
x, w = gauss(2*k + 1)
x = R*x[k:]

# get uniform y grid inside pore
m = int(math.ceil(2.*l0/h))
y = np.linspace(-l0, l0, m)

# list of (x,y) values in pore
Xpore = list(product(x, y))

# uniform x grid outside
l = int(math.ceil(Rx/hout))
x = np.linspace(0., Rx, l)

# gauss y grid outside
L = (Ry-l0-r)
n = int(math.ceil(L/hout))
y, w = gauss(2*n)
yp, ym = l0+r + L*y[n:], -l0-r + L*y[:n]

Xtop = list(product(x, yp))
Xbot = list(product(x, ym))

# TODO: more points at edge y=l0 in x=0,r0

print k+1, m
print l, n
print len(Xpore), len(Xtop), len(Xbot)
