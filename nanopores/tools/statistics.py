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

"""TODO: provide generic/particulare fitting methods like "cdf",
"log_cdf", "moments", "maximum_likelihood", ..."""

_empty = dict()
    
class RandomVariable(object):
    
    parameters = []

    def __init__(self, **params):
        self.params = params
        self.constants = {}
        self.inputs = {}
        self.fixed = dict.fromkeys(self.parameters, False)
        
        for name, X in params.items():
            if not name in self.parameters:
                continue
            if isinstance(X, RandomVariable):
                self.inputs[name] = X
                if X.is_fixed:
                    self.fixed[name] = True
            else:
                self.constants[name] = X
                if X is not None:
                    self.fixed[name] = True
                    
        self.is_fixed = all(self.fixed.values())
        
    def params_sample(self, shape=None):
        params = {}
        for Xname, X in self.inputs.items():
            params[Xname] = X.sample(shape)
        for xname, x in self.constants.items():
            params[xname] = x
        return params
        
    # assuming determined parameters
    def sample(self, shape=None):
        "generate sample of length N"
        params = self.params_sample(shape)
        return self.sample_(shape, **params)
    
    def pdf(self, x, N=100):
        "probability density computed by taking means over input samples"
        shape = [1]*x.ndim + [N]
        params = self.params_sample(shape)
        print params.keys()
        X = self.pdf_(x[..., None], **params)
        return np.mean(X, axis=-1)
    
    # the following should be overloaded to specify model:
    def sample_(self, shape, **params):
        return np.random.random(shape)
    
    def pdf_(self, x, **params):
        return 1.*((0 < x) & (x < 1))
    
    def __repr__(self):
        name = type(self).__name__
        params = ", ".join(["%s=%s" % item for item in self.params.items()])
        return "%s(%s)" % (name, params)
    
    def recursive_params(self):
        params = dict(self.constants)
        for X in self.inputs.values():
            params.update(X.recursive_params())
        return params
    
    def print_params(self):
        print ", ".join([
                "%s=%s" % item for item in self.recursive_params().items()])
    
    def cdf_empirical(self, x, N=100, log=False):
        "sample empirical cdf"
        
    def fit(self, sample, method="cdf", **fit_params):
        fit_function = getattr(self, "fit_" + method)
        fit_function(sample, **fit_params)
        
    def fit_cdf(self, empirical=True, log=False, N=100):
        if empirical:
            cdf = lambda x: self.cdf_empirical(x, log=log, N=N)
            
            
class Poisson(RandomVariable):
    parameters = ["a"]
    
    def sample_(self, shape, a):
        return np.random.poisson(a, shape)
    
class PoissonG0(RandomVariable):
    parameters = ["a"]
    
    # TODO: sample until all are greater zero
    def sample_(self, shape, a):
        k = np.random.poisson(a, shape)
        return k[k > 0]
    
    
class Gamma(RandomVariable):
    parameters = ["K", "tau"]
    
    def sample_(self, shape, K, tau):
        return np.random.gamma(K, scale=tau, size=shape)
    
    def pdf_(self, t, K, tau):
        print K
        return stats.gamma.pdf(t, K, scale=tau)
    
        
if __name__ == "__main__":
    K0 = PoissonG0(a=5) # Providing a value means fixing the parameter
    T0 = Gamma(K=K0, tau=2.) # An RV as parameter generates a compound distr.
    sample = T0.sample(100)
    
    K = PoissonG0(a=None) # None means it has to be fitted from samples
    T = Gamma(K=K, tau=None) # fits a, tau and so determines K and T
    T0.print_params()
    
    t = np.linspace(0., 50., 100)
    plt.hist(sample, normed=True, bins=30)
    plt.plot(t, T0.pdf(t))