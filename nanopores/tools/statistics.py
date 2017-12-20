# (c) 2017 Gregor Mitscha-Baude
"""Generic framework to define, fit, sample, and plot statistical models;
especially suited for complicated compound models, models who include complicated
machinery, and generally models that are not analytically tractable.
Uses simulated annealing as generic parameter search (=fitting) procedure.

Based on the general framework we implement complicated distributions such as
the compound Gamma-Poisson-Gamma-Poisson distribution.

Should enable the following usage:

> K = Poisson(a=None) # None means it has to be fitted from samples
> T = Gamma(K=K, tau=2.0) # Providing a value means fixing the parameter
> T.fit(samples) # fits a and therefore determines both K and T
> k = K.sample(N=1000)
> plt.hist(k, bins=20)
> plt.plot(k, K.pdf(k))
> plt.figure()
> t = np.linspace(0., 5., 100)
> plt.hist(samples)
> plt.plot(t, T.cdf(t))
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

"""TODO: provide generic/particular fitting methods like "cdf",
"log_cdf", "moments", "maximum_likelihood", ..."""

_empty = dict()
    
class RandomVariable(object):
    
    i = 0
    parameters = {}

    def __init__(self, **params):
        self.i = RandomVariable.i
        RandomVariable.i += 1
        self.constants = dict(self.parameters)
        self.inputs = {}
        self.fixed = dict.fromkeys(self.parameters.keys(), False)
        
        for name, X in params.items():
            if not name in self.parameters:
                continue
            if isinstance(X, RandomVariable):
                self.inputs[name] = X
                self.constants.pop(name)
                if X.is_fixed:
                    self.fixed[name] = True
            else:
                if X is not None:
                    self.constants[name] = X
                    self.fixed[name] = True
                    
        self.is_fixed = all(self.fixed.values())
        self.population = self.missing()
        self.shape_pop = ()
    
    def params(self):
        return dict(self.constants, **self.inputs)
    
    # assuming determined parameters
    def sample_params(self, shape):
        params = {}
        for Xname, X in self.inputs.items():
            params[Xname] = X.sample(shape)
        for xname, x in self.constants.items():
            params[xname] = x
        return params
    
    def sample(self, shape=None):
        "generate sample of length N"
        params = self.sample_params(shape)
        return self.sample_(shape, **params)
    
    def pdf(self, x, N=1000, compute_std=True, log=False):
        "probability density computed by taking means over input samples"
        shape = [1]*x.ndim + [N]
        params = self.sample_params(shape)
        X = self.pdf_(x[..., None], **params)
        if compute_std:
            self._std = np.std(X, axis=-1)/np.sqrt(N)
        m = np.mean(X, axis=-1)
        if log:
            m = x * m
            self._std = x * self._std
        return m
    
    def cdf(self, x, N=1000, compute_std=True):
        "cdf computed by taking means over input samples"
        shape = [1]*x.ndim + [N]
        params = self.sample_params(shape)
        X = self.cdf_(x[..., None], **params)
        if compute_std:
            self._std = np.std(X, axis=-1)/np.sqrt(N)
        return np.mean(X, axis=-1)
    
    # versions of the same functions for many parameters at once
    # accessing self.population for missing constants
    # try to be general with x, population shape
    def sample_params_vec(self, shape):
        # assume mean will be over last dimension
        shape_param = shape[:-1] + (1,)
        params = {}
        for Xname, X in self.inputs.items():
            params[Xname] = X.sample_vec(shape)
        for xname, x in self.constants.items():
            if not xname in self.population:
                params[xname] = x
            else:
                params[xname] = np.reshape(self.population[xname], shape_param)
        return params
    
    def sample_vec(self, shape=None):
        params = self.sample_params_vec(shape)
        return self.sample_(shape, **params)

    def pdf_vec(self, x, N=1000, log=False):
        shape_sample = (1,)*x.ndim + self.shape_pop + (N,)
        shape_x = x.shape + (1,)*(len(self.shape_pop) + 1)
        params = self.sample_params_vec(shape_sample)
        X = self.pdf_(x.reshape(shape_x), **params)
        m = np.mean(X, axis=-1)
        if log:
            m = x.reshape(x.shape + (1,)*len(self.shape_pop)) * m
        return m
    
    def cdf_vec(self, x, N=1000):
        shape_sample = (1,)*x.ndim + self.shape_pop + (N,)
        shape_x = x.shape + (1,)*(len(self.shape_pop) + 1)
        params = self.sample_params_vec(shape_sample)
        X = self.cdf_(x.reshape(shape_x), **params)
        return np.mean(X, axis=-1)
    
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
    
    def cdf_(self, x, **params):
        # TODO: make this work
        sample = self.sample_((1000,), **params)
        return empirical_cdf(x, sample)
    
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
    
    def fit(self, sample, method="cdf", **fit_params):
        fit_function = getattr(self, "fit_" + method)
        return fit_function(sample, **fit_params)
        
    def fit_cdf(self, sample, log=False, N=100, **anneal_params):
        "optimize MSE on cdf with annealing."
        xi = grid(sample, tail=0.005, log=log, N=50)
        yi = empirical_cdf(xi, sample)[:, None]
        def F(xi, yi):
            return np.mean((self.cdf_vec(xi, N=N) - yi)**2, 0)
        return self.anneal(F, xi, yi, **anneal_params)
    
    def fit_pdf(self, sample, log=False, N=100, **anneal_params):
        "optimize MSE on cdf with annealing."
        xi = grid(sample, tail=0.005, log=log, N=50)
        xi, yi = empirical_pdf(xi, sample, log=log)
        yi = yi[:, None]
        def F(xi, yi):
            return np.mean((self.pdf_vec(xi, N=N, log=log) - yi)**2, 0)
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
        print "minimum", min(f)
        return self.recursive_missing()
    
    def fit_naive(self, sample):
        params = self.fit_(sample)
        update = {k: params[k] for k in params if k not in self.fixed[k]}
        self.update(update)
        for k in self.inputs:
            if k in update:
                sample = np.array([update[k]])
                self.inputs[k].fit_naive(sample)
    
    def plot_cdf(self, x, *args, **kwargs):
        fx = self.cdf(x, N=1000, compute_std=True)
        itv = 2*self._std
        line, = plt.plot(x, fx, *args, **kwargs)
        plt.fill_between(x, fx - itv, fx + itv,
                         color=line.get_color(), alpha=0.2)
        
    def plot_pdf(self, x, *args, **kwargs):
        log = False if not "log" in kwargs else kwargs.pop("log")
        # plot at centers for compatibility with hist
        x = .5*(x[1:] + x[:-1])
        fx = self.pdf(x, N=100, compute_std=True, log=log)
        itv = 2*self._std
        line, = plt.plot(x, fx, *args, **kwargs)
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
    #return np.array([p[max(0, i-1) : min(n, i+1)].mean() for i in range(n+1)])
    
def smooth(a, k=3):
    a = np.array(a)
    # range of kernel
    start = -(k // 2)
    end = (k // 2) + (k % 2)
    N = a.shape[0]
    b = [a[max(0, i + start) : min(N, i + end)].mean() for i in range(N)]
    return np.array(b)

def grid(data, N=100, tail=0.01, log=False):
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
    "more efficient generation compared to Poisson with condition"
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
#        # solve a0 = a/(1-exp(a)) for a
#        return dict(a=a)
    
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
    
class ScalarMultiple(RandomVariable):
    "implement t*X where X is a RandomVariable and t>0"
    parameters = dict(t=1., X=1.)
        
    def sample_(self, shape, t, X):
        return t*X
    
    def __repr__(self):
        return "%.2g*%s" % (self.constants["t"], repr(self.inputs["X"]))
    
    # TODO: currently can not fit ScalarMultiples because there is no cdf/pdf
#        
#    def sample(self, shape=None):
#        "generate sample of length N"
#        return self.t*self.X.sample(shape)
#        
#    def pdf(self, x, N=1000, compute_std=True):
#        "probability density computed by taking means over input samples"
#        t = self.t
#        mean = self.X.pdf(x/t, N=N, compute_std=compute_std)/t
#        if compute_std:
#            self._std = self.X._std
#        return mean
#    
#    def cdf(self, x, N=1000, compute_std=True):
#        "probability density computed by taking means over input samples"
#        t = self.t
#        mean = self.X.cdf(x/t, N=N, compute_std=compute_std)
#        if compute_std:
#            self._std = self.X._std
#        return mean
#    
#    def sample_vec(self, shape=None):
#        return self.t*self.X.sample_vec(shape)
#    
#    def cdf_vec(self, x, N=1000):
#        t = self.t
#        return self.X.cdf_vec(x/t, N=N)

    
if __name__ == "__main__":
    example1 = True
    
    if example1: # example with stacked Gamma-Poisson-Distributions
        # construct target distribution
        K = ZeroTruncatedPoisson(a=20.) # Assigning a constant means fixing the parameter
        Ta = Gamma(K=K, tau=.005) # An RV as parameter generates a compound distr.
        N = ZeroTruncatedPoisson(a=1e-3*Ta) # RVs can be multiplied by scalar
        T = Gamma(K=N, tau=1e-3)
        # get samples
        sample1 = Ta.sample(1000)
        sample = T.sample(1000)
        
        # construct fitting distribution
        K = ZeroTruncatedPoisson(a=None) # None means it will be fitted
        Ta = Gamma(K=K, tau=None)
        N = ZeroTruncatedPoisson(a=None*Ta)
        T = Gamma(K=N, tau=None)
        
        # fitting methods
        method = "pdf" # or "cdf"; pdf seems to be more robust
        log = True # True should generally yield better fits
        
        # first fit Ta and fix parameters, then fit T
        Ta.fit(sample1, method=method, log=log)
        Ta.fix()
        T.fit(sample, method=method, log=log)
        
        # plot fitted cdf vs. empirical cdf from sample
        plt.figure("cdf")
        Ta.compare_cdfs(sample1, log=True)
        T.compare_cdfs(sample, log=True)
        plt.figure("pdf")
        Ta.compare_pdfs(sample1, log=True)
        T.compare_pdfs(sample, log=True)