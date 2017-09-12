# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from scipy.special import digamma
from scipy.stats import poisson, gamma
from matplotlib import pyplot as plt

euler = 0.577215664901532

t = np.linspace(0, 10, 1000)
plt.plot(t, t*np.exp(t)/np.expm1(t))
plt.show()
exit()
#plt.plot(t, digamma(t))
#plt.plot(t, np.log(t/(1 - np.exp(-t))), ".-")
#


def summand(k, b):
    return digamma(k)*poisson.pmf(k, b)

def f(b, N=50):
    k = np.arange(1, N)
    return np.sum(summand(k, b))*1./np.expm1(b) - np.log(b*np.exp(b)/np.expm1(b))

def finv(x):
    return x/euler + 1

#plt.plot(t, [(f(x) + euler - x*euler + 0.74694*x**2 - 0.336*x**3) for x in t], ".-")
plt.plot(t, [finv(f(x)) for x in t], ".-")
plt.plot(t, t)

#plt.figure()
k = np.arange(1, 20)
#plt.plot(k, summand(k, 0.01), "s--")
plt.show()