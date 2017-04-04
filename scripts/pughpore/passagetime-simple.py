# -*- coding: utf-8 -*-
from __future__ import unicode_literals

# (c) 2017 Gregor Mitscha-Baude
# TODO: obtain rD from actual simulation
from nanopores import fields, kT, eta, qq, savefigs
from numpy import exp, pi, sqrt, linspace, diff, array, dot

L = 46e-9 # length of pore

r = 2.0779e-9 # radius of protein trypsin
V = 0.08 # applied potential
E = V/L # electric field
rD = 0.2
#D = rD* kT/(6.*pi*eta*r) # diffusion constant (Stokes)

# load translocation events without binding
name = "events3_nobind_new"
fields.set_dir_dropbox()
data = fields.get_fields(name)

# take only events that translocated
data.pop("b1")
data.pop("b2")
data, _ = fields._subset(data, data["ood"], lambda x: x==0)
data, times = fields._sorted(data, data["t"])

print "mean"
D = array(data["Dzavg"]).mean() * 1e-9
F = -array(data["Fzavg"]).mean()
v = D/kT * F # electrophoretic velocity
print "D = ", D, "F = ", F, "v = ", v

print "at x = (0,0,0)"
D = 6.8e-12
F = 1.5e-11
v = D/kT * F # electrophoretic velocity
print "D = ", D, "F = ", F, "v = ", v

def mean(lst):
    return sum(lst)/float(len(lst))

def maximum_likelihood(times, n=10):
    times = 1e-3 * array(times)

    T = mean(times)
    Tinv = mean([1./t for t in times])

    def amean(v):
        return mean([1./(1. + L/(v*t)) for t in times])

    def fix(v):
        a = amean(v)
        factor = (sqrt((a-.5)**2  + T*Tinv*a*(1-a)) - (a-.5))/(1-a)
        print a
        #print factor
        return L/T * factor

    v = L*sqrt(Tinv/T) # this initial guess is accurate to 1e-7!!
    for i in range(n):
        v0 = v
        #print "i = %d: v = %s" % (i, v)
        v = fix(v)
        print "i = %d: dv = %s" % (i, abs(v-v0))


    D = v**2/2.*T - v*L + L**2/2.*Tinv
    return v, D

v, D = maximum_likelihood(times)
print "maximum likelihood"
print "D = ", D, "F = ", v*kT/D, "v = ", v

# simple 1D model from Talaga2009
def p(t, timescale=1.):
    # timescale: 1 -> s, 1e-3 -> ms etc
    t *= timescale
    return exp(-(L - t*v)**2/(4.*t*D)) * (L + t*v) / (4.*t * sqrt(pi*t*D))

def pp(times, timescale=1.):
    return array([p(t, timescale) for t in times])

def integrate(t, pt):
    pt = array(pt)
    dt = diff(t)
    values = 0.5*(pt[:-1] + pt[1:])
    return dot(values, dt)

def integrate_hist(hist):
    n, bins, _ = hist
    dt = diff(bins)
    return dot(n, dt)


# scale times
scale = 1e-6 # microseconds
times = [t*1e-3/scale for t in times]

from matplotlib import pyplot as plt
t = linspace(1e-9/scale, 8e-6/scale, 500)
hist = plt.hist(times, bins=30, color="#aaaaff", linewidth=0.5,
                weights=[1./500.]*len(times),
                label="BD simulations")

pt = pp(t, scale) * integrate_hist(hist) * scale
plt.plot(t, pt, "-", color="g", linewidth=3, label="FPT model")
plt.legend(loc="upper right", frameon=False)
plt.xlabel(u"dwell time [µs]")
plt.ylabel(u"rel. frequency")

print "integral", integrate_hist(hist), "==", integrate(t, pt)

#plt.figure()
#plt.hist(data["Fzavg"], bins=30, color="#aaaaff", linewidth=0.5)
#
#plt.figure()
#plt.hist(data["Dzavg"], bins=30, color="#aaaaff", linewidth=0.5)

from folders import FIGDIR
savefigs("current-nobind-hist", FIGDIR + "/rw")

