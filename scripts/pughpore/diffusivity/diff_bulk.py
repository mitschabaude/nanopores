import numpy as np
import os
from math import sinh, acosh
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
import nanopores.tools.fields as fields
from matplotlib import pyplot as plt

geop = nano.Params(pughpore.params)
fields.set_dir(os.path.expanduser("~") + "/Dropbox/nanopores/fields")

r = geop.rMolecule
eps = 1e-8
x = np.linspace(r+eps, 30, 100)

def Cp1(l):
    return 1.-(9./16.)*r/l

def Cp(h):
    x = r/h
    return 1. - 9./16.*x + 1./8.*x**3 - 45./256.*x**4 - 1/16.*x**5

def Cn(l):
    alpha = acosh(l/r)
    s = 0.
    for n in range(1, 100):
        n = float(n)
        K = n*(n+1)/(2*n-1)/(2*n+3)
        s += K*((2*sinh((2*n+1)*alpha)+(2*n+1)*sinh(2*alpha))/(4*(sinh((n+.5)*alpha))**2-(2*n+1)**2*(sinh(alpha))**2) - 1)
    return 1./((4./3.)*sinh(alpha)*s)

Dn = np.array([Cn(xx) for xx in x])
Dt = Cp(x)

def matrix(d):
    return [[d[0], 0., 0.], [0., d[1], 0.], [0., 0., d[2]]]

data = dict(x=list(x), D=list(map(matrix, list(zip(Dn, Dt, Dt)))))
if not fields.exists("pugh_diff_bulk", rMolecule=r):
    print("SAVING...")
    fields.save_fields("pugh_diff_bulk", dict(rMolecule=r), **data)
    fields.update()

plt.plot(x, Dt, ".-", label=r"$D_t$")
plt.plot(x, Cp1(x), ".-", label=r"$D_t$ simple")
plt.plot(x, Dn, ".-", label=r"$D_t$")
plt.legend(loc="lower right")
#plt.show()
