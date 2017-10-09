# (c) 2017 Gregor Mitscha-Baude
from matplotlib import pyplot as plt
import numpy as np

import nanopores
import nanopores.models.randomwalk as rw
from nanopores.tools import fields
fields.set_dir_dropbox()

name = "rw_wei_3"
#data = rw.load_results(name)
#
#rw.histogram(data, a=-2, b=6)
#rw.hist_poisson(data, "attempts", (1, 10))

csvfile = "tau_off_wei.csv"
data = np.genfromtxt(csvfile, delimiter=',')
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

# TODO: need better experimental data => webtool
# now let's reproduce the plot
# first create fake data samples that reproduce the histogram
fake = np.array([])

frac = 1.
while int(counts[0]*frac) > 1:
    frac /= 2.
    a, b = bins[1]*frac, bins[1]*2*frac
    sample = a*(b/a)**(np.random.rand(int(counts[0]*frac)))
    fake = np.append(fake, sample)
    print "frac", frac

for i in range(1, N):
    a, b = bins[i], bins[i+1]
    sample = a*(b/a)**(np.random.rand(counts[i]))
    fake = np.append(fake, sample)

# now get the same number of samples with binding from our own data
data = rw.load_results(name)
bind = data.bindings > 0
times = 1e-9*data.times # data are in ns, we want s

n = sum(counts)
print "Found %d simulated binding events, have %d experimental binding events." % (sum(bind), n)
if sum(bind) > n:
    # only consider same number of bind_times
    i = np.flatnonzero(bind)[n-1] + 1
else:
    i = len(times)
print "Number of simulated events we use:", i
bind_times = times[:i][bind[:i]]
fail_times = times[:i][data.fail[:i]]
success_times = times[:i][data.success[:i]]

# simple histogram
plt.figure("hist_simple")
plt.hist(fake, bins=bins, label="Wei et al. 2012")
plt.hist(bind_times, bins=bins, histtype="step", rwidth=1., label="Simulation")
plt.legend()
plt.xlabel(r"$\tau$ off [s]")
plt.ylabel("Count")
plt.xlim(0, 20)

plt.figure("hist_attempts")
rw.hist_poisson(data, "attempts", (0, 12), lines=False)
plt.figure("hist_bindings")
rw.hist_poisson(data, "bindings", (0, 4), lines=False)
#rw.hist_poisson(data, "attempts", (1, 10), modified=True)

# histogram plot with short and long events and experimental data
plt.figure("hist")
cutoff = 0.03e-3 # cutoff frequency in s
a, b = -6.5, 3 # log10 of plot interval
bins = np.logspace(a, b, 40)

# successful events
hist = plt.hist(success_times, bins=bins, color="green", rwidth=0.9, label="Translocated")
# failed attempts
hist = plt.hist(fail_times, bins=bins, color="red", rwidth=0.9, label="Did not translocate")


#total = rw.integrate_hist(hist, cutoff)
#tmean = times[times > cutoff].mean()
#T = np.logspace(a-3, b, 1000)
#fT = np.exp(-T/tmean)*T/tmean
#fT *= total/integrate_values(T, fT, cutoff)
#plt.plot(T, fT, label="exp. fit, mean = %.2f ms" % (tmean,),
#         color="dark" + color, **params)
#plt.xlim(10**a, 10**b)

#fake0 = fake[fake > bins[1]]
plt.hist(fake, bins=bins, histtype="step", color="orange", label="Wei et al. 2012")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Count")
plt.xlabel(r"$\tau$ off [s]")
plt.ylim(ymin=1., ymax=1e5)
plt.legend()
#plt.show()

import folders
nanopores.savefigs("tau_off", folders.FIGDIR + "/wei", (4, 3))