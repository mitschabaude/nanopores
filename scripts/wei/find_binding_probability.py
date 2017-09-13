# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from matplotlib import pyplot as plt
import nanopores
import nanopores.models.randomwalk as randomwalk
from nanopores.tools import fields
from nanopores.tools.solvers import cache_forcefield
fields.set_dir_dropbox()

params = nanopores.user_params(
    # general params
    geoname = "wei",
    dim = 2,
    rMolecule = 1.25, # 6.
    h = 5.,
    Nmax = 1e5,
    Qmol = 2., #15.,
    bV = -0.2,
    dp = 23.,

    # random walk params
    N = 1000, # number of (simultaneous) random walks
    dt = 1., # time step [ns]
    walldist = 2., # in multiples of radius, should be >= 1
    margtop = 50.,
    margbot = 0.,
    zstart = 46.5, # 46.5
    xstart = 0., # 42.
    rstart = None,

    # receptor params
    rec_t = 3.7e9,
    rec_p = 0.01763,
    rec_eps = 0.1,
)

receptor = randomwalk.Ball([8.5 - 3. + 2., 0., 40.5], 0.5) # ztop 46.5
receptor_params = dict(
    exclusion = False,
    walldist = 1.,
    #minsize = 0.01, # accuracy when performing reflection

    binding = True,
    eps = params.rec_eps, # margin in addition to walldist, determines re-attempting [nm]
    t = params.rec_t, # mean of exponentially distributed binding duration [ns]
    p = params.rec_p, # binding probability for one attempt

    use_force = False, # if True, t_mean = t*exp(-|F|*dx/kT)
    dx = 0.1, # width of bond energy barrier [nm]
)

name = "rw_wei_p"

@cache_forcefield(name, default=params)
def binding_prob(P, **params):
    for p, result in nanopores.collect_dict(P):
        pore = nanopores.get_pore(**params)
        rw = randomwalk.RandomWalk(pore, **params)
        receptorp = dict(receptor_params)
        receptorp["p"] = p
        rw.add_domain(receptor, **receptorp)

        print "Start Random Walk with p = %s." % p
        for t in rw.walk():
            pass

        result.new = dict(
            p0 = np.mean(rw.bindings > 0),
            a0 = np.mean(rw.attempts > 0),
            mean_bindings = np.mean(rw.bindings),
            mean_attempts = np.mean(rw.attempts),
            std_bindings = np.std(rw.bindings),
            std_attempts = np.std(rw.attempts),
            mean_time = np.mean(rw.times),
            std_time = np.std(rw.times),
            mean_log_time = np.mean(np.log(rw.times)),
        )
    return result

def binding_prob_from_data(rMolecule=1.25, rPore=6.):
    phys = nanopores.Physics()

    # calculate binding probability with data from (Wei 2012)
    kon = 20.9e6 # association rate constant [1/Ms] = binding events per second
    c = 180e-9 # concentration [M = mol/l = 1000 mol/m**3]
    cmol = c * 1e3 * phys.mol # concentration [1/m**3]
    ckon = c*kon

    # Smoluchowski rate equation gives number of arrivals at pore entrance per sec
    D = phys.kT / (6. * phys.pi * phys.eta * rMolecule * 1e-9) # [m**2/s]
    r = rPore* 1e-9 # effective radius for proteins at pore entrance [m]
    karr = 2.*phys.pi * r * D * cmol # arrival rate
    p0 = ckon / karr # fraction of events that have binding

    print "Average time between events (tau_on): %.2f s (from experimental data)" % (1./ckon)
    print "Number of bindings per second: %.1f (inverse of tau_on)" % ckon
    print "Number of events per second: %.1f (from Smoluchowski rate equation)" % karr
    print "=> fraction of events where binding occurs: %.1f / %.1f = %.5f" % (ckon, karr, p0)
    print "= p0 = prob. of binding at least once"
    return p0

def invert_monotone(y, X, Y):
    """assuming data pairs X, Y = f(X) coming from an increasing function,
    find x such that f(x) = y, by interpolating between the nearest y values."""
    # get index of first element in Y larger than y
    i = np.searchsorted(Y, y)
    # now interpolate linearly between (X[i-1], Y[i-1]) and (X[i], Y[i])
    # y = Y[i-1] + t*(Y[i] - Y[i-1])
    t = (y - Y[i-1])/(Y[i] - Y[i-1])
    x = X[i-1] + t*(X[i] - X[i-1])
    return x

if __name__ == "__main__":
    #P1 = np.linspace(0, 1, 50)
    #data1 = binding_prob(P1, nproc=5, N=1000)
    #plt.plot(P1, data1.p0, ".-", label="actual (N=1000)")

    P2 = np.linspace(0, 1, 100)
    data2 = binding_prob(P2, nproc=5, N=20000)
    plt.plot(P2, data2.p0, ".-", label="actual (N=20000)")

    P3 = np.linspace(0, 0.05, 10)
    data3 = binding_prob(P3, nproc=5, N=100000)
    #plt.plot(P3, data3.p0, ".-", label="actual (N=100000)")

    PP = np.linspace(0, 1, 500)
    plt.plot(PP, 1. - np.exp(-0.3*PP), label="Poisson")

    p0 = binding_prob_from_data()
    plt.plot([0, 1], [p0, p0], "k--", label="p0 = %.4f" % p0)

    p = invert_monotone(p0, P2, data2.p0)
    plt.plot([p], [p0], "ok", label="inferred p = %.4f" % p)

    print "binding prob. inferred from simulations: p = %.6f" % p
    a = 0.3
    ap = -np.log(1 - p0)
    p1 = ap/a
    print "binding prob. inferred from assumed Poisson distribution: p = %.6f" % p1

    plt.xlabel("p")
    plt.ylabel("probability of >= 1 binding")
    plt.legend()
    plt.show()