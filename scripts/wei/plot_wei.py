# (c) 2017 Gregor Mitscha-Baude
from matplotlib import pyplot as plt
import numpy as np

import nanopores
import nanopores.models.randomwalk as rw
from nanopores.tools import fields
fields.set_dir_dropbox()

name = "rw_wei_2"
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

# now let's reproduce the plot
# first create fake data samples that reproduce the histogram
fake = np.array([])
for i in range(N):
    a, b = bins[i], bins[i+1]
    sample = a + (b-a)*np.random.rand(counts[i])
    fake = np.append(fake, sample)

plt.hist(fake, bins=bins, alpha=0.3, label="Wei et al. 2012")

# now get the same number of samples with binding from our own data
data = rw.load_results(name)
bind = data.bindings > 0
times = data.times[bind]
times *= 1e-9 # data are in ns, we want s

n = sum(counts)
print "Found %d simulated binding events, have %d experimental binding events." % (sum(bind), n)
if sum(bind) > n:
    times = times[:n]

plt.hist(times, bins=bins, alpha=0.3, label="Simulation")
plt.legend()
plt.xlabel(r"$\tau$ off [s]")
plt.ylabel("count")
plt.xlim(0, 20)

rw.hist_poisson(data, "attempts", (1, 10))
rw.hist_poisson(data, "bindings", (1, 10))

plt.show()