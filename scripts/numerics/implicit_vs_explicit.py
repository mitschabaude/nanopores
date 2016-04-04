"compare implicit and explicit molecule model visually for 2D Howorka pore"

import os, numpy, Howorka
from dolfin import *
from nanopores import *
#from matplotlib.pyplot import figure, plot, legend, show, title, xlabel, ylabel, savefig

add_params(
himp = .2,
hexp = .2,
Nimp = 1e4,
Nexp = 1e4,
z0 = 2.,
Qmol = -1.,
)

# simulate explicit molecule
def explicit():
    geo, phys = Howorka.setup2D(z0=z0, h=hexp, Qmol=Qmol)
    #plot(geo.pwconst("permittivity"))
    pb, pnps = Howorka.solve2D(geo, phys, Nmax=Nexp, cheapest=True)
    (v, cp, cm, u, p) = pnps.solutions()
    pnps.visualize("fluid")
    return geo, phys, v, u
        
# simulate implicit molecule
def implicit():
    geo, phys = Howorka.setup2D(z0=None, h=himp)
    pb, pnps = Howorka.solve2D(geo, phys, Nmax=Nimp, cheapest=True)
    (v, cp, cm, u, p) = pnps.solutions()
    #pnps.visualize("fluid")
    return v, u
              
geo, phys, vexp, uexp = explicit()
vimp, uimp = implicit()
    
# analytical potential from molecule in homogeneous dielectric fluid
eps = phys.eperm*phys.rpermw
q = phys.qq
lscale = phys.lscale
x0 = [0., z0]

class vAna(Expression):
    def eval(self, value, x):
        r = sqrt(sum((t - t0)**2 for t, t0 in zip(x, x0)))
        value[0] = -q*lscale/(4.*pi*eps*r)   

synonyme = dict(
    notmol = {"fluid", "dna", "membrane"}
)
geo.import_synonymes(synonyme)
mesh = geo.submesh("notmol")
V = FunctionSpace(mesh, "CG", 1)

v = Function(V)
va = Function(V)
v0 = Function(V)
dv = Function(V)

v.interpolate(vexp)
v0.interpolate(vimp)
plot(v, title="explicit")
plot(v0, title="implicit")

va.interpolate(vAna())
plot(va, title="analytical")

dv.vector()[:] = v.vector()[:] - v0.vector()[:] #- va.vector()[:]
plot(dv, title="difference")

U = VectorFunctionSpace(mesh, "CG", 1)

u = Function(U)
u0 = Function(U)
du = Function(U)

u.interpolate(uexp)
u0.interpolate(uimp)
plot(u, title="explicit")
plot(u0, title="implicit")

du.vector()[:] = u.vector()[:] - u0.vector()[:]
plot(du, title="difference")
interactive()
