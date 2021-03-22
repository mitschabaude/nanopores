# (c) 2017 Gregor Mitscha-Baude
from matplotlib import pyplot as plt
import numpy as np

csvfile = "tau_on_wei.csv"
data = np.genfromtxt(csvfile, delimiter=',')
bins = data[:, 0]
counts = data[:, 1]

# inspection showed that there seems to be a good,
# evenly spaced approximation to all bins except the first with
# spacing 0.4, i.e. of the form (beta + 0.4*np.arange(0, N)) for some beta
dx = 0.04
x = bins[:-1]
N = len(x)
# minimize norm(x - (beta + 0.55*np.arange(0, N)) w.r.t. beta
beta = x.mean() - dx*(N-1)/2.
print("beta", beta)
# beta is close to 0.016, which we take
bins = 0.016 + dx*np.arange(0, N)
bins = [0.] + list(bins)

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
plt.xlim(0, 1.5)
plt.show()