# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import nanopores
from find_binding_probability import (binding_prob,
                                      binding_prob_from_data, invert_monotone)

# load data
P2 = np.linspace(0, 1, 100)
P2a = P2[P2 > 0.05]
data2 = binding_prob(P2, nproc=5, calc=False, N=20000)
data2a = binding_prob(P2a, nproc=5, calc=False, N=20000)

P3 = np.linspace(0, 0.05, 10)
data3 = binding_prob(P3, nproc=5, calc=False, N=100000)

P4a = np.linspace(0.01, 0.5, 20)
data4a = binding_prob(P4a, nproc=5, calc=False, N=4000)

P4 = np.linspace(0.01, 0.5, 20)
data4 = binding_prob(P4, nproc=5, calc=True, N=100000)

a = 0.3 # mean number of attempts
a1 = 2.2
p0 = binding_prob_from_data()
p = invert_monotone(p0, P3, data3.p0)
pmod = a/a1
#pmod = p0/(1. - np.exp(-a1*p))

# big plot
plt.figure("p0")
plt.plot(P2, data2.p0, ".", label="Simulated (N=20000)")
PP = np.linspace(0, 1, 500)
plt.plot(PP, 1. - np.exp(-a*PP), label="Poisson")
plt.xlabel(r"$p$")
plt.ylabel(r"$p_0$") #probability of >= 1 binding")
plt.legend(frameon=False)

# smaller plot where p is inferred
plt.figure("p0_small")
#plt.plot(P2, data2.p0, "o", label="Simulated (N=20000)", zorder=100)
plt.plot(P3, data3.p0, "o", label="Simulated (N=100000)", zorder=100)
PP = np.linspace(0, 1, 500)
plt.plot(PP, 1. - np.exp(-0.3*PP), label="Poisson (a = 0.3)")
plt.plot(PP, pmod*(1. - np.exp(-a1*PP)), label="Mod. Poisson (a = 2.2)")
plt.plot([0, 1], [p0, p0], "k--", label="p0 from data")

plt.plot([p], [p0], "o", color="#000000", label="inferred p = %.3f" % p, zorder=100)
#plt.axvline(x=p, ymin=0., ymax=p0/0.025, color="#000000", zorder=-90)

plt.xlim(-0.002, 0.062)
plt.ylim(-0.001, 0.023)
plt.yticks([0, .005, .01, .015, .02])
plt.xlabel(r"Binding probability $p$")
plt.ylabel(r"$p_0$") #probability of >= 1 binding")
plt.legend(frameon=False)

# big plot
plt.figure("p0_fit")
plt.plot(P2, data2.p0, ".", label="Simulated (N=20000)", zorder=100)
PP = np.linspace(0, 1, 500)
plt.plot(PP, 1. - np.exp(-a*PP), label="Poisson (a = 0.3)")
plt.plot(PP, pmod*(1. - np.exp(-a1*PP)), label="Mod. Poisson (a = 2.2)")
print "pmod", pmod
plt.xlabel(r"Binding probability $p$")
plt.ylabel(r"$p_0$") #probability of >= 1 binding")
plt.gca().add_patch(Rectangle((-0.01, -0.002), 0.07, 0.02, fc="none", ec="k"))
plt.legend(frameon=False)

import folders
nanopores.savefigs("binding_prob", folders.FIGDIR + "/wei", (4, 3))

print "binding prob. inferred from simulations: p = %.6f" % p
ap = -np.log(1 - p0)
p1 = ap/a
print "binding prob. inferred from assumed Poisson distribution: p = %.6f" % p1

# plot mean, std, log of time
plt.figure("time_stdmean")
mu = np.array(data4.mean_time)
sigma = np.array(data4.std_time)
log = np.array(data4.mean_log_time)
plt.plot(P4, (mu/sigma)**2, "o", label="Simulated (N=100000)")
def f(x):
    return np.exp(x)/(2.*np.expm1(x)/x - 1.)
plt.plot(P4, f(a*P4), label="Poisson (a = %.1f)" % a)
plt.plot(P4, f(a1*P4), label="Mod. Poisson (a = %.1f)" % a1)
plt.legend(frameon=False)

plt.figure("time_log")
euler = 0.577215664901532
theta = -0.573810187498/euler + 1. # estimate from histogram
plt.plot(P4, (log - np.log(mu))/euler + 1., "o", label="Simulated (N=100000)")
plt.plot(P4, np.ones_like(P4)*theta, "--k", label="Estimate from histogram")
from scipy.special import digamma
from scipy.stats import poisson

def summand(k, b):
    return digamma(k)*poisson.pmf(k, b)
def f1(b, N=50):
    k = np.arange(1, N)
    return np.sum(summand(k, b))*1./np.expm1(b) - np.log(b*np.exp(b)/np.expm1(b))
def f1v(x):
    return np.array([f1(b) for b in x])

plt.plot(P4, f1v(a*P4)/euler + 1., label="Poisson (a = %.1f)" % a)
plt.plot(P4, f1v(a1*P4)/euler + 1., label="Mod. Poisson (a = %.1f)" % a1)
plt.xlabel("p")
plt.legend(frameon=False)

plt.figure("time_mean")
#tau = 3.7
tau = 3.88
def g(x):
    return x/(1. - np.exp(-x))
def taufit(a):
    return tau/g(a*p)

plt.plot(P4, 1e-9*taufit(a1)/tau*mu, "o", label="Simulated (N=100000)", zorder=100)
plt.plot(P4, tau*np.ones_like(P4), "--", color="orange", label="Const., tau = %.2f" % tau)
plt.plot(P4, taufit(a)*g(a*P4), label="Poisson, tau = %.2f" % (taufit(a)), color="C1")
plt.plot(P4, taufit(a1)*g(a1*P4), label="Mod. Poisson, tau = %.2f" % (taufit(a1)), color="C2")
plt.plot([p], [tau], "o", color="#000066", label=r"p, tau off from data")
#lima, limb = plt.ylim()
#plt.axvline(x=p, ymin=0., ymax=(tau - lima)/(limb - lima), color="#000066", zorder=-90)
#plt.xlim(-0.01, 0.21)
plt.ylim(ymax=6.6)
plt.xlabel(r"Binding probability $p$")
plt.ylabel(r"Mean $\tau$ off [s]")
plt.legend(loc="upper left", frameon=False)

plt.figure("time_std")
sig = 4.02
def h(x):
    return np.sqrt(g(x)*(2. - x/np.expm1(x)))
def sigfit(a):
    return sig/h(a*p)

plt.plot(P4, 1e-9*sigfit(a1)/tau*sigma, "o", label=r"Simulated (N=100000)")
plt.plot(P4, sigfit(a)*h(a*P4), label=r"Poisson (a = %.1f, $\tau$ = %.2f)" % (a, sigfit(a)))
plt.plot(P4, sigfit(a1)*h(a1*P4), label=r"Mod. Poisson (a = %.1f, $\tau$ = %.2f)" % (a1, sigfit(a1)))
plt.plot(P4, sig*np.ones_like(P4), "--", color="orange", label=r"Exponential ($\tau$ = %.2f)" % sig)

plt.plot([p], [sig], "o", color="#000066", label="p, sigma inferred from data")
#lima, limb = plt.ylim()
#plt.axvline(x=p, ymin=0., ymax=(tau - lima)/(limb - lima), color="#000066", zorder=-90)
plt.xlabel("p")
plt.ylabel("std. dev. of binding duration [s]")
plt.legend(frameon=False)

import folders
nanopores.savefigs("tau", folders.FIGDIR + "/wei", (5, 3.7))
#plt.show()