# (c) 2017 Gregor Mitscha-Baude
# TODO: obtain rD from actual simulation
from nanopores import fields, kT, eta, qq, savefigs
from numpy import exp, pi, sqrt, linspace, diff, array, dot

L = 46e-9 # length of pore

r = 2.0779e-9 # radius of protein trypsin
v = 0.08 # applied potential
E = v/L # electric field

rD = 0.2
D = rD* kT/(6.*pi*eta*r) # diffusion constant (Stokes)
v = D/kT * 5*qq*E # electrophoretic velocity

# simple 1D model from Talaga2009
def p(t):
    return exp(-(L - t*v)**2/(4.*t*D)) * (L + t*v) / (4.*t * sqrt(pi*t*D))
    
def pp(times):
    return array([p(t) for t in times])
    
def integrate(t, pt):
    pt = array(pt)
    dt = diff(t)
    values = 0.5*(pt[:-1] + pt[1:])
    return dot(values, dt)
    
def integrate_hist(hist):
    n, bins, _ = hist
    dt = diff(bins)
    return dot(n, dt)
    
# load translocation events without binding
name = "events_nobind"
fields.set_dir_dropbox()
data = fields.get_fields(name)

times = sorted(data["t"])

from matplotlib import pyplot as plt
t = linspace(1e-8, 250e-6, 500)
hist = plt.hist(times, bins=30, color="#aaaaff", linewidth=0.5, label="BD simulations")

pt = pp(t) * integrate_hist(hist)
plt.plot(t, pt, "-", color="g", linewidth=3, label="FPT model")
plt.legend(loc="upper right")

from folders import FIGDIR
savefigs("current-nobind-hist", FIGDIR + "/rw")

