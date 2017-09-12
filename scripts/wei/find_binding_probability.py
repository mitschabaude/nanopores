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
            mean_bindings = np.mean(rw.bindings),
            mean_attempts = np.mean(rw.attempts),
            std_bindings = np.std(rw.bindings),
            std_attempts = np.std(rw.attempts),
        )
    return result

if __name__ == "__main__":
    P = np.linspace(0, 1, 50)
    data = binding_prob(P, nproc=5, N=1000)
    plt.plot(P, data.p0, ".-", label="actual (N=1000)")

    P = np.linspace(0, 1, 100)
    data = binding_prob(P, nproc=5, N=20000)
    plt.plot(P, data.p0, ".-", label="actual (N=20000)")

    PP = np.linspace(0, 1, 500)
    plt.plot(PP, 1. - np.exp(-0.3*PP), label="Poisson")
    plt.xlabel("p")
    plt.ylabel("probability of >= 1 binding")
    plt.legend()
    plt.show()