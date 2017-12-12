# (c) 2017 Gregor Mitscha-Baude
"""Random walks in Pugh pore on 2D proxy geometry, to determine distribution of
attempt time needed to fit binding parameters."""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import tangent
import dolfin
import nanopores
import nanopores.models.randomwalk as randomwalk
from nanopores.tools.polygons import Rectangle
from nanopores.tools import fields
fields.set_dir_mega()

params = nanopores.user_params(
    # general params
    # geo
    geoname = "pughcyl",
    dim = 2,
    diamPore = 6.,
    rMolecule = 2.0779,
    R = 40.,
    Htop = 60.,
    Hbot = 35.,
    geop = dict(R=40., Htop=60., Hbot=35.),
    x0 = None,
    # physics
    Qmol = 5.,
    bulkcon = 1000.,
    dnaqsdamp = 0.7353,
    bV = -0.1,
    posDTarget = True,
    # solver
    h = 2.,
    frac = 0.5,
    Nmax = 5e4,
    imax = 30,
    tol = 1e-3,
    cheapest = False,
    stokesiter = False, #True
    hybrid = True,
    reconstruct = False,

    # random walk params
    N = 100, # number of (simultaneous) random walks
    dt = .2, # time step [ns]
    walldist = 1., # in multiples of radius, should be >= 1
    margtop = 15., # 35 is about the maximum
    margbot = 0.,
    rstart = 5.,
    initial = "sphere",

    # receptor params
    bind_everywhere = False,
    lbind = 8., # length of binding zone for long binding [nm]
    ra = .2, # binding zone radius [nm] (1nm means whole channel)
    collect_stats_mode = True,
)
########### WHAT TO DO  ###########
todo = nanopores.user_params(
    test_solver = False,
    plot_dolfin = False,
    plot_streamlines = False,
    video = False,
    plot_distribution = True,
)

########### SETUP ###########
NAME = "rw_pugh_0"

def binding_(params):
    # TODO: make other params dependent on additional **args
    return dict(
    binding = True,
    bind_type = "zone",
    collect_stats_mode = params["collect_stats_mode"],
    t = 1e9, # mean of exponentially distributed binding duration [ns]
    ka = 1e9, # (bulk) association rate constant [1/Ms]
    ra = params["ra"], # radius of the association zone [nm]
    use_force = True, # if True, t_mean = t*exp(-|F|*dx/kT)
    dx = 0.1, # width of bond energy barrier [nm]
    )

def setup_rw(params, **_):
    pore = nanopores.get_pore(**params)
    rw = randomwalk.RandomWalk(pore, **params)
    binding_params = binding_(params)
    if params["bind_everywhere"]:
        rw.add_wall_binding(**binding_params)
    else:
        r = pore.params.l3/2.
        z = -pore.params.hpore/2. + pore.params.h4
        wbind = 1.
        lbind = params["lbind"]
        bindsite = Rectangle((r, r + wbind), (z, z + lbind))
        rw.add_domain(bindsite, exclusion=False, **binding_params)
        
    #rw.domains[0].__dict__.update(exclusion=True, binding=True,
    #          bind_type="collision", eps=0., p=1., walldist=1.)
    return rw

########### TEST AND PLOT SOLVER ###########
if todo.test_solver:
    import nanopores.models.nanopore as nanopore
    setup = nanopore.Setup(**params)
    _, pnps = nanopore.solve(setup, True)    
    dolfin.interactive()

if todo.plot_dolfin:
    rw = setup_rw(params)
    dolfin.plot(rw.D[1])
    dolfin.plot(rw.F)
    dolfin.interactive()
    
if todo.plot_streamlines:
    rw = setup_rw(params)
    rw.plot_streamlines(both=True, R=20, Hbot=30, Htop=35,
                        maxvalue=1e-10, figsize=(5, 5))
    plt.figure("D")
    dolfin.plot(rw.D[1], backend="matplotlib")
    #plt.show()
    
if todo.video:
    rw = setup_rw(params)
    randomwalk.run(rw)

########### PLOT ATTEMPT TIME DISTRIBUTION ###########
def NLS(ti, yi, t0=0., tol=1e-14):
    "nonlinear least squares to fit exponential distribution with mean exp(x)"
    xi = np.log(ti)
    
    def f(x, xi, yi):
        return np.sum((1. - np.exp(-np.exp(xi - x)) - yi)**2)
    
    # minimize f by solving df(x) = 0 with newton method
    df = tangent.grad(f)
    ddf = tangent.grad(df)
    x = np.log(t0)
    while np.abs(df(x, xi, yi)) > tol:
        x -= df(x, xi, yi)/ddf(x, xi, yi)
        #print "|f(x)|", np.abs(df(x, xi, yi))
    return np.exp(x)

def NLS2(ti, yi, t10=1., t20=100., w0=0.5, tol=1e-14):
    "nonlinear least squares to fit DOUBLE exponential distribution"
    xi = np.log(ti)
    # find: characteristic times exp(x1), exp(x2) and weights w1, w2
    
    def f(theta, xi, yi):
        w = theta[0]
        x1, x2 = theta[1], theta[2]
        z = w*np.exp(-np.exp(xi - x1)) + (1.-w)*np.exp(-np.exp(xi - x2))
        return np.sum((1. - z - yi)**2)
    
    # minimize f by solving df(x) = 0 with newton method
    df = tangent.grad(f)
    ddf = tangent.grad(df, mode="forward")
    def Jf(theta, xi, yi):
        return np.array([ddf(theta, xi, yi, 1., [1,0,0]),
                         ddf(theta, xi, yi, 1., [0,1,0]),
                         ddf(theta, xi, yi, 1., [0,0,1])])
    
    theta = np.array([w0, np.log(t10), np.log(t20)])
    dftheta = df(theta, xi, yi)
    while np.linalg.norm(dftheta) > tol:
        print "|(grad f)(theta)|", np.linalg.norm(dftheta)
        theta -= np.linalg.solve(Jf(theta, xi, yi), dftheta)
        dftheta = df(theta, xi, yi)
        
    return theta[0], np.exp(theta[1]), np.exp(theta[2])

def NLS_general(F, xi, yi, p0=1., tol=1e-12):
    "nonlinear least squares to fit arbitrary f with any number of parameters"
    # F = F(x, p), p0 MUST match len(p), F takes array of x    
    # F must be compatible with tangent module
    def f(p, xi, yi):
        return np.sum((F(xi, p) - yi)**2)
    
    # minimize f by solving df(x) = 0 with newton method
    n = len(p0)
    df = tangent.grad(f)
    ddf = tangent.grad(df, mode="forward")
    ei = lambda i: np.eye(1, n, i)[0, :]
    
    def Jf(p, xi, yi):
        return np.array([ddf(p, xi, yi, 1., ei(i)) for i in range(n)])
    
    p = np.array(p0)
    dfp = df(p, xi, yi)
    while np.linalg.norm(dfp) > tol:
        print "|grad f|", np.linalg.norm(dfp)
        p -= np.linalg.solve(Jf(p, xi, yi), dfp)
        dfp = df(p, xi, yi)
        
    return tuple(p)

def NLS_bruteforce(F, xi, yi, p, width=1., N=100):
    # make p 1d array
    p = np.atleast_1d(p)
    # create parameter range
    # TODO: this currently only applies to len(p)==2
    x = np.logspace(-width, width, N)
    xx = np.column_stack((np.repeat(x[:, None], len(x), 0),
                          np.tile(x[:, None], [len(x), 1])))
    pp = p[None, :] * xx
    
    f = np.sum((F(xi[None, :], pp) - yi)**2, 1)
    i = np.argmin(f)
    print "minimum:", f[i]
    print "parameters:", pp[i, :]
    return tuple(pp[i, :])

def NLS_annealing(F, xi, yi, p, sigma=5., N=100, n=10):
    # sigma = initial (multiplicative) standard deviation
    # N = size of population in one iteration
    # n = number of iterations
    print "initial", p
    p = np.atleast_1d(p)
    dim = len(p)
    factor = 0.5 # to reduce sigma per iteration
    # make initial sigma act like multiplication by sigma^(+-1)
    sigma = np.log(sigma)*np.ones(dim)
    
    for k in range(n):
        # create new population by adding multiplicative gaussian noise
        P = p[None, :] * np.exp(np.random.randn(N, dim) * sigma[None, :])
        # compute mean square loss on population
        f = np.mean((F(xi[None, :], P) - yi)**2, 1)
        # replace p by new best guess
        p = P[np.argmin(f), :]
        # update sigma
        sigma *= factor
        print "parameters:", p
    print "minimum", min(f)
        
    return tuple(p)
        
def fit_gamma(ti):
    mu = ti.mean()
    sigma = ti.std()
    C = mu**2/sigma**2
    
    def f(x):
        return np.exp(x)/(2.*np.expm1(x)/x - 1.)

    def f1(x):
        y = 2.*np.expm1(x)/x - 1.
        z = 2./x**2*(x*np.exp(x) - np.expm1(x))
        return np.exp(x)*(1 - z/y)/y
    
    # newton solve
    K = 2.*C # initial value
    #print "Newton iteration:"
    for i in range(10):
        dK = -(f(K) - C)/f1(K)
        K = K + dK
    #    print i, "Residual", f(K) - C, "Value K =", K
    #print
    
    tau = mu*(1. - np.exp(-K))/K
    return K, tau

from scipy.special import iv
class CompoundGamma(object):
    
    def __init__(self, ti):
        self.ti = ti
        self.K, self.tau = self.fit_brute(ti)
        
    def fit_todo(self, ti, n=40):
        pass
        
    def fit_naive(self, ti):
        return fit_gamma(ti)
    
    def fit_brute(self, ti, n=100):
        # first guess to get orders of magnitude right
        p = np.array(fit_gamma(ti))
        # define problem # TODO:
        #xi = np.sort(ti)[int(len(ti)/n/2)::int(len(ti)/n)]
        #yi = np.arange(len(xi))/float(len(xi))
        bins = np.logspace(np.log10(min(ti)), np.log10(max(ti)), n)
        hist, _ = np.histogram(ti, bins=bins)
        xi = 0.5*(bins[:-1] + bins[1:])
        yi = np.cumsum(hist)/float(np.sum(hist))
        # minimize
        #K, tau = NLS_bruteforce(self.cfd_vec, xi, yi, p, width=1., N=100)
        K, tau = NLS_annealing(self.cfd_vec, xi, yi, p, N=100, n=20)
        return K, tau
    
    def pdf_direct(self, tt, N=50):
        a = self.K
        t = tt/self.tau
        S = np.ones_like(t)
        s = np.ones_like(t)
        for k in range(1, N):
            s *= (a*t)/(k*(k+1.))
            S += s
        return np.exp(-t)*a/np.expm1(a) * S /self.tau
    
    def pdf_bessel(self, tt):
        a = self.K
        t = tt/self.tau
        return np.exp(-t)*np.sqrt(a/t)*iv(1., 2.*np.sqrt(a*t))/self.tau/np.expm1(a)    
    
    def cdf(self, tt, N=50):
        a = self.K
        tau = self.tau
        gamma = sp.stats.gamma.cdf
        S = np.zeros_like(tt)
        s = 1.
        for k in range(1, N):
            s *= a/k
            S += s*gamma(tt, k, scale=tau)
        return 1./np.expm1(a) * S
    
    def cfd_vec(self, tt, p, N=50):
        # cdf that takes parameter as vector input, for fitting
        a = p[:, 0:1]
        tau = p[:, 1:2]
        gamma = sp.stats.gamma.cdf
        S = np.ones((p.shape[0], tt.shape[1]))
        s = np.ones((1, tt.shape[0]))
        for k in range(1, N):
            s = s*a/k
            S = S + s*gamma(tt, k, scale=tau)
        return 1./np.expm1(a) * S
    
    def gammainc(self, tt, k,tau):
        # TODO: can not differentiate wrt k
        # implement fitting yourself
        #for j in range(k):
        pass
        

if todo.plot_distribution:
    N = 10000
    params0 = dict(params, N=N)
    rw = randomwalk.get_rw(NAME, params0, setup=setup_rw)
    ta = rw.attempt_times
    ta1 = ta[ta > 0.]
    tmean = ta1.mean()
    
    bins = np.logspace(np.log10(min(ta1)), np.log10(max(ta1)), 35)
    #bins = np.logspace(-3., 2., 35)
    hist, _ = np.histogram(ta1, bins=bins)
    cfd = np.cumsum(hist)/float(np.sum(hist))
    t = 0.5*(bins[:-1] + bins[1:])
    
    #n = 100
    #t = np.sort(ta1)[int(len(ta1)/n/2)::int(len(ta1)/n)]
    #cfd = np.arange(len(t))/float(len(t))
    
    plt.figure("ta_cfd", figsize=(4,3))
    tt = np.logspace(np.log10(min(ta1)), np.log10(max(ta1)), 100)
    plt.semilogx(t, cfd, "v", label="Simulations")
    
    # naive exp. fit
    #plt.semilogx(tt, 1. - np.exp(-tt/tmean), label="Simple exp. fit")
    
    # proper exp. fit
    toff = NLS(t, cfd, t0=tmean)
    #plt.semilogx(tt, 1. - np.exp(-tt/toff), label="Exp. fit")
    
    # ML gamma fit
    #K, _, tau = sp.stats.gamma.fit(ta1)
    K = (tmean/ta1.std())**2
    tau = tmean/K
    #plt.semilogx(tt, sp.stats.gamma.cdf(tt, K, scale=tau), label="Simple Gamma fit")
    gamma = CompoundGamma(ta1)
    plt.semilogx(tt, gamma.cdf(tt), label="Compound Gamma fit")
    
    # double exponential fit
    #w, toff1, toff2 = NLS2(t, cfd, t10=toff/2., t20=toff*2., w0=.4)
    #plt.semilogx(tt, 1. - w*np.exp(-tt/toff1) - (1.-w)*np.exp(-tt/toff2),
    #             label="Double exp. fit")
    
    plt.xlabel("Attempt time [ns]")
    plt.ylabel("Cumulative frequency")
    plt.xlim(xmin=1.)
    plt.legend()
    
    xlog = False
    plt.figure("ta_hist", figsize=(4,3))
    #tt = np.logspace(-1, 4., 20)
    bins = np.linspace(0., 200., 30)
    tt = np.linspace(0., 200., 300)
    plt.figure("ta_hist", figsize=(4,3))
    plt.hist(ta1, bins=bins, normed=True, log=False, label="Simulations")
    
    #tt0 = tt if xlog else 1.
    #plt.plot(tt, tt0/tmean * np.exp(-tt/tmean),
    #         label="Simple exp. fit, mean=%.3gns" % tmean)
    #plt.plot(tt, tt0/toff * np.exp(-tt/toff),
    #         label="Exp. fit, mean=%.3gns" % toff)
    #dt = rw.dt
    #kk = np.arange(1000)
    #k0 = tmean/dt
    
    #plt.plot(tt, sp.stats.gamma.pdf(tt, K, scale=tau),
    #         label="Simple Gamma fit")
    plt.plot(tt, gamma.pdf_direct(tt), label="Compound Gamma fit")
    #plt.plot(tt, gamma.pdf_bessel(tt), ".k", label="Compound Gamma fit")
    #plt.plot(kk*dt, poisson.pmf(kk, k0), label="Poisson fit")
    #plt.plot(tt)
    if xlog:
        plt.xscale("log")
    #plt.ylim(ymin=1e-10)
    #plt.yscale("log")
    plt.xlabel("Attempt time [ns]")
    plt.ylabel("Rel. frequency")
    plt.legend()
    

import folders
nanopores.savefigs("rw_cyl", folders.FIGDIR + "/pugh", ending=".pdf")