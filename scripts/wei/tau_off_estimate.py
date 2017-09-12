# (c) 2017 Gregor Mitscha-Baude
import numpy as np

# get experimental event data
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

# estimate b and tau parameters
euler = 0.577215664901532
theta = np.mean(np.log(fake)) - np.log(np.mean(fake))
print "theta:", theta
b = theta/euler + 1.
print "estimate b:", b
b0 = b/np.expm1(b)*np.exp(b)
print "mean number of bindings, given K>0:", b0
tau = np.mean(fake)/b0
print "mean binding duration, exponential assumption:", np.mean(fake)
print "mean binding duration, gamma assumption:", tau
