# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from matplotlib import pyplot as plt
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

P4 = np.linspace(0.01, 0.5, 20)
data4 = binding_prob(P4, nproc=5, calc=False, N=4000)

a = 0.3 # mean number of attempts
a1 = 2.17

# big plot
plt.figure("p0")
plt.plot(P2, data2.p0, ".-", label="actual (N=20000)")
PP = np.linspace(0, 1, 500)
plt.plot(PP, 1. - np.exp(-a*PP), label="Poisson")
plt.xlabel(r"$p$")
plt.ylabel(r"$p_0$") #probability of >= 1 binding")
plt.legend()

# smaller plot where p is inferred
plt.figure("p0_small")
plt.plot(P2a, data2a.p0, ".-", label="actual (N=20000)")
plt.plot(P3, data3.p0, ".-", color="darkblue", label="actual (N=100000)")
PP = np.linspace(0, 1, 500)
plt.plot(PP, 1. - np.exp(-0.3*PP), label="Poisson")

p0 = binding_prob_from_data()
plt.plot([0, 1], [p0, p0], "k--", label="p0 = %.4f" % p0)

p = invert_monotone(p0, P3, data3.p0)
plt.plot([p], [p0], "o", color="#000066", label="inferred p = %.4f" % p)
plt.axvline(x=p, ymin=0., ymax=p0/0.025, color="#000066", zorder=-90)

plt.xlim(0, 0.08)
plt.ylim(0, 0.025)
plt.xlabel(r"$p$")
plt.ylabel(r"$p_0$") #probability of >= 1 binding")
plt.legend()

# big plot
plt.figure("p0_fit")
plt.plot(P2, data2.p0, ".-", label="actual (N=20000)")
PP = np.linspace(0, 1, 500)
plt.plot(PP, 1. - np.exp(-a*PP), label="Poisson (a = 0.3)")
plt.plot(PP, p0*(1. - np.exp(-a1*PP))/(1. - np.exp(-a1*p)), label="Mod. Poisson (a = 2.2)")
pmod = p0/(1. - np.exp(-a1*p))
print "pmod", pmod
plt.xlabel(r"$p$")
plt.ylabel(r"$p_0$") #probability of >= 1 binding")
plt.legend()

print "binding prob. inferred from simulations: p = %.6f" % p
a = 0.3
ap = -np.log(1 - p0)
p1 = ap/a
print "binding prob. inferred from assumed Poisson distribution: p = %.6f" % p1

# plot mean, std, log of time
plt.figure("time_stdmean")
mu = np.array(data4.mean_time)
sigma = np.array(data4.std_time)
log = np.array(data4.mean_log_time)
plt.plot(P4, (mu/sigma)**2, ".-", label="actual (N=4000)")
def f(x):
    return np.exp(x)/(2.*np.expm1(x)/x - 1.)
plt.plot(PP[1:], f(a*PP[1:]), label="Poisson")
plt.legend()

plt.figure("time_log")
euler = 0.577215664901532
theta = -0.736486800755/euler + 1. # estimate from histogram
plt.plot(P4, (log - np.log(mu))/euler + 1., ".-", label="actual (N=4000)")
plt.plot(P4, np.ones_like(P4)*theta, "--k", label="estimate from histogram")
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
plt.plot(P4, f1v(a1*P4)/euler + 1., label="Poisson (a = %.1f)" % a1)
plt.legend(loc="lower right")

plt.figure("time_mean")
tau = 3.7
def g(x):
    return x/(1. - np.exp(-x))
def taufit(a):
    return tau/g(a*p)

a1 = 2.2
plt.plot(P4, 1e-9*taufit(a1)/tau*mu, ".-", label="actual (N=4000)")
plt.plot(P4, taufit(a)*g(a*P4), label="Poisson (a = %.1f, tau = %.2f)" % (a, taufit(a)))
plt.plot(P4, taufit(a1)*g(a1*P4), label="Poisson (a = %.1f, tau = %.2f)" % (a1, taufit(a1)))
plt.plot(P4, tau*np.ones_like(P4), "--", color="orange", label="Exponential (tau = %.2f)" % tau)

plt.plot([p], [tau], "o", color="#000066", label="p, tau inferred from data")
#lima, limb = plt.ylim()
#plt.axvline(x=p, ymin=0., ymax=(tau - lima)/(limb - lima), color="#000066", zorder=-90)
plt.ylabel("p")
plt.ylabel("mean binding duration [s]")
plt.legend()

plt.figure("time_std")
sig = 4.02
def h(x):
    return np.sqrt(g(x)*(2. - x/np.expm1(x)))
def sigfit(a):
    return sig/h(a*p)

a1 = 2.2
plt.plot(P4, 1e-9*sigfit(a1)/tau*sigma, ".-", label="actual (N=4000)")
plt.plot(P4, sigfit(a)*h(a*P4), label="Poisson (a = %.1f, tau = %.2f)" % (a, sigfit(a)))
plt.plot(P4, sigfit(a1)*h(a1*P4), label="Poisson (a = %.1f, tau = %.2f)" % (a1, sigfit(a1)))
plt.plot(P4, sig*np.ones_like(P4), "--", color="orange", label="Exponential (tau = %.2f)" % sig)

plt.plot([p], [sig], "o", color="#000066", label="p, sigma inferred from data")
#lima, limb = plt.ylim()
#plt.axvline(x=p, ymin=0., ymax=(tau - lima)/(limb - lima), color="#000066", zorder=-90)
plt.ylabel("p")
plt.ylabel("std. dev. of binding duration [s]")
plt.legend()

import folders
nanopores.savefigs("binding_prob", folders.FIGDIR + "/wei", (6, 4.5))
#plt.show()