# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from numpy import exp, expm1
from matplotlib import pyplot as plt

def f(x):
    return exp(x)/(2.*expm1(x)/x - 1.)

def f1(x):
    y = 2.*expm1(x)/x - 1.
    z = 2./x**2*(x*exp(x) - expm1(x))
    return exp(x)*(1 - z/y)/y

def solve(C, n=20):
    x = 2.*C # initial value
    print "Newton iteration:"
    for i in range(n):
        dx = -(f(x) - C)/f1(x)
        x = x + dx
        print i, "Residual", f(x) - C, "Value", x
    print
    return x

def get_parameters(mu, sigma):
    C = mu**2/sigma**2
    ap = solve(C, 10)
    lmbda = ap/(mu*(1. - exp(-ap)))
    return ap, lmbda

# get experimental event data
csvfile = "tau_off_wei.csv"
data = np.genfromtxt(csvfile, delimiter=",")
bins = data[:, 0]
counts = data[:, 1]

# inspection showed that there seems to be a good,
# evenly spaced approximation to all bins except the first and last with
# spacing 0.55, i.e. of the form (beta + 0.55*np.arange(0, N)) for some beta
x = bins[:-1]
N = len(x)
# minimize norm(x - (beta + 0.55*np.arange(0, N)) w.r.t. beta
beta = x.mean() - 0.55*(N-1)/2.
# turns out beta is close to 0.25, which gives nice numbers,
# so we will just take that
bins = 0.25 + 0.55*np.arange(0, N)
bins = [0.] + list(bins) + [20.]
N = N+1

# the counts should be integer-values, so
counts = np.round(counts).astype(int)

# now let's reproduce the plot
# first create fake data samples that reproduce the histogram
fake = np.array([])

frac = 1.
while frac > 0.5: #int(counts[0]*frac) > 1:
    frac /= 2.
    a, b = bins[1]*frac, bins[1]*2*frac
    sample = a*(b/a)**(np.random.rand(int(counts[0]*frac)))
    fake = np.append(fake, sample)
    print "frac", frac

for i in range(1, N):
    a, b = bins[i], bins[i+1]
    sample = a*(b/a)**(np.random.rand(counts[i]))
    fake = np.append(fake, sample)

# compute mu, variance, solve for parameters
mu = np.mean(fake)
sigma = np.std(fake)
ap, lmbda = get_parameters(mu, sigma)

print "mu, sigma =", mu, sigma
print "mu^2/sigma^2 =", mu**2/sigma**2
print "ap, lambda =", ap, lmbda
print
print "binding probability: %.3f (for a=2.2)" % (ap/2.2,)
print "mean binding duration: %.3f s" % (1./lmbda,)

#
#
#X = linspace(0., 10., 100)
#plt.plot(X, f(X, 0.))
#plt.plot(X, X/2.)
#plt.show()
