"create 2D point set for Howorka model where force shall be evaluated."
import numpy as np
from itertools import product
import math
import matplotlib.pyplot as plt
import nanopores
gauss = np.polynomial.legendre.leggauss

nanopores.add_params(
    h = 0.5,
    hout = 1.,    
    Ry = 10.,
    Rx = 3.,
)

def points(h, hout, r0, r, l0, Rx, Ry):
    # effective pore radius for molecule
    R = r0 - r
    
    # get scaled 1D gauss quadrature nodes
    # for x grid inside pore
    # this has h = h/2 to resolve the thin region
    k = int(math.ceil(R/h*2.))
    x, w = gauss(2*k + 1)
    x = R*x[k:]
    
    # get uniform y grid inside pore
    m = int(math.ceil(2.*l0/h))
    y = np.linspace(-l0-r+h/2, l0+r-h/2, m)
    
    # list of (x,y) values in pore
    Xpore = list(product(x, y))
    
    # gauss x grid outside
    l = int(math.ceil(Rx/h))
    x, w = gauss(2*l)
    x = Rx-Rx*x[l:]
    #x = np.linspace(0., Rx, l)
    
    # gauss y grid outside
    L = (Ry-l0-r)
    n = int(math.ceil(L/hout))
    y, w = gauss(2*n)
    yp = l0+r + L + L*y[:n]
    ym = -yp
    
    Xtop = list(product(x, yp))
    Xbot = list(product(x, ym))
    
    # combine
    X = Xtop + Xpore + Xbot
    return X
    # TODO: more points at edge y=l0 in x=0,r0 ?

# load parameters and create points
from nanopores.tools import fields
from nanopores.geometries.H_cyl_geo.params_geo import r0, rMolecule, l0, r1

X = points(h, hout, r0, rMolecule, l0/2, r1 + rMolecule, Ry)
x, y = [z[0] for z in X], [z[1] for z in X]
plt.scatter(x, y)
plt.show()

#fields.save_entries("xforce", PARAMS, X=X, N=len(X))
#print "Created and stored %d evaluation points." %(len(X),)
#fields.update()