# (c) 2017 Gregor Mitscha-Baude
"""Generic framework to define, fit, sample, and plot statistical models;
especially suited for complicated compound models, models who include complicated
machinery, and generally models that are not analytically tractable.
Uses simulated annealing as generic parameter search (=fitting) procedure.

Based on the general framework we implement complicated distributions such as
the compound Gamma-Poisson-Gamma-Poisson distribution.

Example usage:

> K = Poisson(a=None) # None means it has to be fitted from samples
> T = Gamma(K=K, tau=2.0) # Providing a value means fixing the parameter
> T.fit(samples) # fits a and therefore determines both K and T
>
> t = np.linspace(0., 5., 100)
> plt.hist(samples, normed=True)
> plt.plot(t, T.pdf(t))
"""

# TODOs:
# -) naive fitting to arrive at better initial guess
# -) check what is missing for discrete RVs (i.e. fitting)
#
# Further ideas:
# -) direct integration for more exact dfs when no. parameters is small
# -) any other stuff leveraging analytical knowledge, e.g. moments
# -) more advanced evolution search strategies that are more reliable
# -) think about what is the common generalization of this and neural nets with
#    backpropagation. strategies for backpropagation combined with sampling?

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

class RandomVariable(object):
    
    i = 0
    parameters = {}
    derived_from = []

    def __init__(self, **params):
        self.i = RandomVariable.i
        RandomVariable.i += 1
        self.constants = dict(self.parameters)
        self.inputs = {}
        self.fixed = dict.fromkeys(self.parameters.keys(), True)
        self.is_derived_from = {k: (True if k in self.derived_from else False
                                 ) for k in self.parameters}
        
        for name, X in params.items():
            if not name in self.parameters:
                continue
            if isinstance(X, RandomVariable):
                self.inputs[name] = X
                self.constants.pop(name)
                if not X.is_fixed:
                    self.fixed[name] = False
            else:
                if X is not None:
                    self.constants[name] = X
                else:
                    self.fixed[name] = False
                    
        self.is_fixed = all(self.fixed.values())
        self.population = self.missing()
        self.shape_pop = ()
        
    def __getattr__(self, attr):
        return self.params()[attr]
    
    def params(self):
        return dict(self.constants, **self.inputs)

    def sample_params(self, shape, train=False):
        params = {}            
        for name, X in self.inputs.items():
            if self.is_derived_from[name]:
                params[name] = X.sample_params(shape, train)
            else:
                params[name] = X.sample(shape, train)
        for name, x in self.constants.items():
            if train and name in self.population:
                params[name] = self.population[name].reshape(shape[:-1] + (1,))
            else:
                params[name] = x
        return params
    
    def sample(self, shape=None, train=False):
        "generate sample of length N"
        params = self.sample_params(shape, train)
        return self.sample_(shape, **params)

    def pdf(self, x, N=1000, train=False, log=False, compute_std=False):
        "probability density computed by taking means over input samples"
        shape_p = self.shape_pop if train else ()
        shape = (1,)*x.ndim + shape_p + (N,)
        x = x.reshape(x.shape + (1,)*(len(shape_p) + 1))
        params = self.sample_params(shape, train)
        
        X = self.pdf_(x, **params)
        factor = x[..., 0] if log else 1. 
        if compute_std:
            self._std = factor * np.std(X, axis=-1)/np.sqrt(N)
        return factor * np.mean(X, axis=-1)
    
    def cdf(self, x, N=1000, train=False, compute_std=False):
        "cdf computed by taking means over input samples"
        shape_p = self.shape_pop if train else ()
        shape = (1,)*x.ndim + shape_p + (N,)
        x = x.reshape(x.shape + (1,)*(len(shape_p) + 1))
        params = self.sample_params(shape, train)
        
        X = self.cdf_(x, **params)
        if compute_std:
            self._std = np.std(X, axis=-1)/np.sqrt(N)
        return np.mean(X, axis=-1)
    
    def fit(self, sample, method="cdf", **fit_params):
        fit_function = getattr(self, "fit_" + method)
        return fit_function(sample, **fit_params)
        
    def fit_cdf(self, sample, log=False, N=100, Ngrid=50, **anneal_params):
        "optimize MSE on cdf with annealing."
        xi = grid(sample, tail=0.005, log=log, N=Ngrid)
        yi = empirical_cdf(xi, sample)[:, None]
        def F(xi, yi):
            fxi = self.cdf(xi, N=N, train=True)
            return np.mean((fxi - yi)**2, 0)
        return self.anneal(F, xi, yi, **anneal_params)
    
    def fit_pdf(self, sample, log=False, N=100, Ngrid=50, **anneal_params):
        "optimize MSE on cdf with annealing."
        xi = grid(sample, tail=0.005, log=log, N=Ngrid)
        xi, yi = empirical_pdf(xi, sample, log=log)
        yi = yi[:, None]
        def F(xi, yi):
            fxi = self.pdf(xi, N=N, train=True, log=log)
            return np.mean((fxi - yi)**2, 0)
        return self.anneal(F, xi, yi, **anneal_params)
        
    def anneal(self, F, xi, yi, n_pop=100, n_it=20, sigma=5., factor=0.5,
               verbose=True):
        "minimize loss of yi = F(xi; p) wrt p with simulated annealing"
        # n_pop = size of population in one iteration
        # n_it = number of iterations
        # sigma = initial (multiplicative) standard deviation
        # factor = factor to reduce sigma per iteration
        if verbose:
            t = self.recursive_missing()
            print "   ".join(map(lambda t: "%s%d" % t[::2], t))
            print " ".join(map(lambda t: "%.2f" % t[1], t))
        for k in range(n_it):
            # create new population by adding multiplicative gaussian noise
            self.spawn_population_lognormal(n_pop=n_pop, sigma=sigma)
            # compute loss
            f = F(xi, yi)
            # replace p by new best guess
            self.update_from_population(np.argmin(f))
            # update sigma
            sigma *= factor
            # print params
            if verbose:
                print " ".join(map(lambda t: "%.2g" % t[1],
                                   self.recursive_missing()))
        if verbose:
            print "minimum", min(f)
        return min(f) #self.recursive_missing()
    
    def fit_naive(self, sample):
        params = self.fit_(sample)
        update = {k: params[k] for k in params if k not in self.fixed[k]}
        self.update(update)
        for k in self.inputs:
            if k in update:
                sample = np.array([update[k]])
                self.inputs[k].fit_naive(sample)
    
    def missing(self):
        return {k: v for k, v in self.constants.items() if not self.fixed[k]}
    
    def fix(self):
        self.fixed = dict.fromkeys(self.fixed.keys(), True)
        self.is_fixed = True
        for X in self.inputs.values():
            X.fix()
    
    def update(self, **params):
        for k in params:
            if k in self.constants:
                self.constants[k] = params[k]
                
    def update_from_population(self, i):
        "set constants to i-th value of population"
        self.update(**{k: v[i] for k, v in self.population.items()})
        for X in self.inputs.values():
            if not X.is_fixed:
                X.update_from_population(i)
                
    def set_population(self, **params):
        first = params.values()[0]
        if np.isscalar(first):
            self.shape_pop = ()
            assert all(np.isscalar(p) for p in params.values())
        else:
            self.shape_pop = first.shape
            assert all(p.shape == self.shape_pop for p in params.values())
        missing = self.missing()
        self.population = {k: params[k] for k in params if k in missing}
                
    def spawn_population_lognormal(self, n_pop=100, sigma=5.):
        # TODO: in principle, uncertain parameters could be RVs themselves
        # then one would have either to discriminate two kinds of inputs
        # or collapse sampling and population spawning
        self.shape_pop = (n_pop,)
        self.population = {}
        for k in self.missing():
            c = self.constants[k]
            self.population[k] = c * np.exp(np.random.randn(n_pop)*sigma)
        for X in self.inputs.values():
            if not X.is_fixed:
                X.spawn_population_lognormal(n_pop, sigma)
    
    # the following 3 should be overloaded to specify model:
    def sample_(self, shape, **params):
        return np.random.random(shape)
    
    # default behaviour is to construct empirical dfs from sample
    #def pdf_(self, t, K, tau):
    #    return stats.gamma.pdf(t, K, scale=tau)
    
    def __repr__(self):
        name = type(self).__name__
        params = self.constants.items() + self.inputs.items()
        params = ", ".join(["%s=%s" % item for item in params])
        return "%s(%s)" % (name, params)
    
    def __mul__(self, t):
        #assert np.isscalar(t)
        return ScalarMultiple(t=t, X=self)
    
    def __rmul__(self, t):
        #assert np.isscalar(t)
        return ScalarMultiple(t=t, X=self)
    
    def recursive_params(self):
        # returns items for nonuniqueness
        params = [(k, v, self.i) for k, v in self.constants.items()]
        for X in self.inputs.values():
            params.extend(X.recursive_params())
        return params
    
    def recursive_missing(self):
        params = [(k, v, self.i) for k, v in self.missing().items()]
        for X in self.inputs.values():
            params.extend(X.recursive_missing())
        return params
    
    def print_params(self):
        print ", ".join([
                "%s=%s (%d)" % item for item in self.recursive_params()])
    
    def plot_cdf(self, x, *args, **kwargs):
        std = True if not "std" in kwargs else kwargs.pop("std")
        fx = self.cdf(x, N=1000, compute_std=std)
        line, = plt.plot(x, fx, *args, **kwargs)
        if std:
            itv = 2*self._std
            plt.fill_between(x, fx - itv, fx + itv,
                         color=line.get_color(), alpha=0.2)
        
    def plot_pdf(self, x, *args, **kwargs):
        log = False if not "log" in kwargs else kwargs.pop("log")
        std = True if not "std" in kwargs else kwargs.pop("std")
        # plot at centers for compatibility with hist
        x = .5*(x[1:] + x[:-1])
        fx = self.pdf(x, N=1000, compute_std=std, log=log)
        line, = plt.plot(x, fx, *args, **kwargs)
        if std:
            itv = 2*self._std
            plt.fill_between(x, fx - itv, fx + itv,
                         color=line.get_color(), alpha=0.2)
        
    def compare_cdfs(self, sample, log=True):
        t = grid(sample, 20, 0.005, log=log)
        tt = grid(sample, 100, 0.005, log=log)
        plt.plot(t, empirical_cdf(t, sample), "o")
        self.plot_cdf(tt)
        if log:
            plt.xscale("log")
            
    def compare_pdfs(self, sample, log=True):
        t = grid(sample, 20, 0.005, log=log)
        tt = grid(sample, 100, 0.005, log=log)
        t, epdf = empirical_pdf(t, sample, log=log)
        plt.plot(t, epdf, "o")
        self.plot_pdf(tt, log=log)
        if log:
            plt.xscale("log")
            
def empirical_cdf(x, data):
    "evaluate edcf of data at points x, i.e. np.mean(data < x) for all x"
    data = np.sort(data)
    return np.searchsorted(data, x)/float(data.size)

def empirical_pdf(x, data, log=False):
    "evaluate epdf at bin centers by creating normed histogram"
    p, _ = np.histogram(data, bins=x)
    mass = np.dot(np.diff(np.log(x)), p) if log else np.dot(np.diff(x), p)
    x = .5*(x[1:] + x[:-1])
    return x, p/mass

def empirical_cdf_vec(x, data):
    # x and data just have to be broadcastable, sampling dimension is last
    return np.mean(data <= x, axis=-1)
    
def smooth(a, k=3):
    a = np.array(a)
    # range of kernel
    start = -(k // 2)
    end = (k // 2) + (k % 2)
    N = a.shape[0]
    b = [a[max(0, i + start) : min(N, i + end)].mean() for i in range(N)]
    return np.array(b)

def grid(data, N=100, tail=0.01, log=False, xmin=None, xmax=None):
    "regularly spaced evaluation nodes spanning data distribution"
    data = np.sort(data)
    n = data.size
    # set limits near tail ends on both sides
    xmin = data[int(np.floor(tail*n))]
    xmax = data[int(np.ceil((1.-tail)*n))-1]
    if log:
        return np.logspace(np.log10(xmin), np.log10(xmax), N)
    else:
        return np.linspace(xmin, xmax, N)

####### commonly used RVs #########
        
class Empirical(RandomVariable):
    
    def __init__(self, data):
        self.data = data
        RandomVariable.__init__(self)
    
    def sample_(self, shape):
        return np.random.choice(self.data, size=shape)
    
    def cdf_(self, x):
        return empirical_cdf(x, self.data)

class Bernoulli(RandomVariable):
    parameters = dict(p=.5)
    
    def sample_(self, shape, p):
        return stats.bernoulli.rvs(p, size=shape)
        
#class Categorical(RandomVariable):
#    """generalization of Bernoulli variable, output is integer 0,...,n-1 with
#    resp. probability p_0,...,p_n-1. n is no parameter, but fixed at
#    instantiation by the length of the probabilities vector."""
#    
#    # TODO: maybe have to initialize .parameters
#    def __init__(self, *p):
#        self.n = len(p)
#        keys = ["p%d" % i for i in range(self.n)]
#        params = dict(zip(keys, p))
#        RandomVariable.__init__(self, **params)
#        
#    def prob(self, pdict):
#        "return sorted and normalized probabilities from dict"
#        p = np.array([x for _, x in sorted(pdict.items())])
#        return p/np.sum(p, axis=0)[None, ...]
#        
#    # FIXME: does not work for vector-valued probabilities
#    def sample_(self, shape, **p):
#        return np.random.choice(self.n, size=shape, p=self.prob(p))
                     
class Poisson(RandomVariable):
    parameters = dict(a=1.)
    
    def sample_(self, shape, a):
        return np.random.poisson(a, shape)
    
    def pdf_(self, x, a):
        return stats.poisson.pmf(x, a)
    
    def cdf_(self, x, a):
        return stats.poisson.cdf(x, a)
    
def broadcast_mask(M):
    "broadcast boolean mask to index possibly larger array"
    colon = slice(None)
    return tuple(i if d > 1 else colon for d, i in zip(M.shape, M.nonzero()))
    
class ZeroTruncatedPoisson(RandomVariable):
    "Poisson conditioned on K > 0"
    parameters = dict(a=1.)
    
    def sample_(self, shape, a):
        if np.isscalar(a) or a.size == 1:
            return self.sample_scalar(shape, a)
        
        #print shape, a.shape
        AMAX = 30
        k = 1
        K = np.full(shape, k)
        
        # for large a values, just use non-truncated poisson
        large = broadcast_mask(a > AMAX)
        small = broadcast_mask(a <= AMAX)
        K[large] = np.random.poisson(a[large], K[large].shape)
        Ksmall = K[small]
        a = a[small]
        
        # actual algorithm begins here
        s = a/np.expm1(a)
        S = s
        U = np.random.random(Ksmall.shape)
        new = S < U
        while np.any(new):
            k += 1
            Ksmall[new] = k
            s = s*a/float(k)
            S = S + s
            new = S < U
        K[small] = Ksmall
        return K
    
    def sample_scalar(self, shape, a):
        AMAX = 30
        if a > AMAX:
            return np.random.poisson(a, shape)
        k = 1
        K = np.full(shape, k)
        s = a/np.expm1(a)
        S = s
        U = np.random.random(shape)
        new = S < U
        while np.any(new):
            k += 1
            K[new] = k
            s = s*a/float(k)
            S = S + s
            new = S < U
        return K
    
    def pdf_(self, x, a):
        return stats.poisson.pmf(x, a)*np.exp(a)/np.expm1(a)
    
    def cdf_(self, x, a):
        return stats.poisson.cdf(x, a)*np.exp(a)/np.expm1(a)

    # TODO:    
#    def fit_(self, data):
#        a0 = np.mean(data)
#        # solve a0 = a/(1-exp(-a)) for a
#        return dict(a=a)

class Exponential(RandomVariable):
    parameters = dict(tau=1.)
    
    def sample_(self, shape, tau):
        return np.random.exponential(scale=tau, size=shape)
    
    def pdf_(self, t, tau):
        return stats.expon.pdf(t, scale=tau)
    
    def cdf_(self, t, tau):
        return stats.expon.cdf(t, scale=tau)
    
    def fit_(self, data):
        return dict(tau=data.mean())  
    
class LeftTruncatedExponential(RandomVariable):
    parameters = dict(tau=1., tmin=0.1)
    
    def sample_(self, shape, tau, tmin):
        umax = np.exp(-tmin/tau)
        u = umax * np.random.random(size=shape)
        return -tau * np.log(u)
    
    def pdf_(self, t, tau, tmin):
        return np.exp(-(t-tmin)/tau)/tau * (1.*(t > tmin))
    
    def cdf_(self, t, tau, tmin):
        return 1. - np.exp(-np.maximum(t - tmin, 0.)/tau)
    
class Gamma(RandomVariable):
    parameters = dict(K=1, tau=1.)
    
    def sample_(self, shape, K, tau):
        return np.random.gamma(K, scale=tau, size=shape)
    
    def pdf_(self, t, K, tau):
        return stats.gamma.pdf(t, K, scale=tau)
    
    def cdf_(self, t, K, tau):
        return stats.gamma.cdf(t, K, scale=tau)
    
    def fit_(self, data):
        K, _, tau = stats.gamma.fit(data)
        return dict(K=K, tau=tau)
    
class LeftTruncatedGamma(RandomVariable):
    parameters = dict(K=1, tau=1., tmin=0.1)
    
    pass
    
####### RVs derived from others #########
        
class ScalarMultiple(RandomVariable):
    "implement t*X where X is a RandomVariable and t > 0"
    parameters = dict(t=1., X=1.)
    derived_from = ["X"]
        
    def sample_(self, shape, t, X):
        x = self.inputs["X"].sample_(shape, **X)
        return t*x
        
    def pdf_(self, x, t, X):
        return self.inputs["X"].pdf_(x/t, **X)/t
    
    def cdf_(self, x, t, X):
        return self.inputs["X"].cdf_(x/t, **X)

    def __repr__(self):
        return "%.2g*%s" % (self.constants["t"], repr(self.inputs["X"]))
    
class OneOf(RandomVariable):
    """RV is one of X, Y where X is w times more likely than Y.
    In other words, [X, Y][i] where i ~ Bernoulli(1/(1+w))"""
    # fitting is biased to X being more likely, to get a unique order
    parameters = dict(X=1., Y=1., w=5.)
    derived_from = ["X", "Y"]
    
    def sample_(self, shape, X, Y, w):
        x = self.X.sample_(shape, **X)
        y = self.Y.sample_(shape, **Y)
        chose_y = np.bool_(stats.bernoulli.rvs(1./(1. + w), size=shape))
        x[chose_y] = y[chose_y]
        return x
    
    def pdf_(self, x, X, Y, w):
        p = 1./(1. + w)
        return (1. - p)*self.X.pdf_(x, **X) + p*self.Y.pdf_(x, **Y)
    
    def cdf_(self, x, X, Y, w):
        p = 1./(1. + w)
        return (1. - p)*self.X.cdf_(x, **X) + p*self.Y.cdf_(x, **Y)
    
def DoubleExponential(tau1=1., tau2=1., w=None):
    #w = p/(1.-p) if p is not None else None
    X = Exponential(tau=tau1)
    Y = Exponential(tau=tau2)
    return OneOf(X=X, Y=Y, w=w)

def LeftTruncatedDoubleExponential(tau1=1., tau2=1., w=None, tmin=0.5):
    #w = p/(1.-p) if p is not None else None
    X = LeftTruncatedExponential(tau=tau1, tmin=tmin)
    Y = LeftTruncatedExponential(tau=tau2, tmin=tmin)
    return OneOf(X=X, Y=Y, w=w)
    
if __name__ == "__main__":
    example1 = True
    example2 = True
    
    if example1: # example with stacked Gamma-Poisson-Distributions
        # construct target distribution
        K = ZeroTruncatedPoisson(a=20.) # Assigning a constant means fixing the parameter
        Ta = Gamma(K=K, tau=0.005) # An RV as parameter generates a compound distr.
        N = ZeroTruncatedPoisson(a=100.*Ta) # RVs can be multiplied by scalar
        T = 1e-3*Gamma(K=N)
        # get samples
        sample1 = Ta.sample(1000)
        sample = T.sample(1000)
        
        # construct fitting distribution
        K = ZeroTruncatedPoisson(a=None) # None means it will be fitted
        Ta = Gamma(K=K, tau=None)
        N = ZeroTruncatedPoisson(a=None*Ta)
        N.a.constants["t"] = 50. # intiial guess t=1. would not work
        T = None*Gamma(K=N)
        
        # fitting methods
        method = "pdf" # or "cdf"; pdf seems to be more robust
        log = True # True should generally yield better fits
        
        # first fit Ta and fix parameters, then fit T
        Ta.fit(sample1, method=method, log=log)
        Ta.fix()
        T.fit(sample, method=method, log=log)
        
        # alternatively, simply treat Ta as coming from a fixed empirical dist.
        Ta_alt = Empirical(sample1)
        N_alt = ZeroTruncatedPoisson(a=None*Ta_alt)
        T_alt = None*Gamma(K=N_alt)
        T_alt.fit(sample, method=method, log=log)
        
        # alternatively, just fit an Exponential
        T_exp = None*Exponential()
        T_exp.fit(sample, method=method, log=log)
        
        # plot fitted cdf vs. empirical cdf from sample
        tt = grid(sample, 100, 0.005, log=True)
        plt.figure("cdf")
        Ta.compare_cdfs(sample1, log=True)
        T.compare_cdfs(sample, log=True)
        T_alt.plot_cdf(tt, ":k")
        T_exp.plot_cdf(tt, "--b")
        plt.figure("pdf")
        Ta.compare_pdfs(sample1, log=True)
        T.compare_pdfs(sample, log=True)
        T_alt.plot_pdf(tt, ":k", log=True)
        T_exp.plot_pdf(tt, "--b", log=True)
        
        print "\nT", T
        print "\nT_alt", T_alt
        print "\nT_exp", T_exp
        
    if example2: # combinations of Exponential variables
        tmin = 0.005
        sample = LeftTruncatedDoubleExponential(
                     tau1=0.01, tau2=1., w=1., tmin=tmin).sample(10000)
        
        T = LeftTruncatedDoubleExponential(
                     tau1=None, tau2=None, w=None, tmin=tmin)
        
        T.fit(sample, method="cdf", log=True, sigma=2., factor=0.9, n_it=50)
        plt.figure("cdf")
        T.compare_cdfs(sample, log=True)
        plt.figure("pdf")
        T.compare_pdfs(sample, log=True)
        
        print T