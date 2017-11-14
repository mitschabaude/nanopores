# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import nanopores
import nanopores.models.randomwalk as randomwalk
from nanopores.tools import fields
fields.set_dir_mega()

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
    N = 100, # number of (simultaneous) random walks
    dt = 1., # time step [ns]
    walldist = 2., # in multiples of radius, should be >= 1
    margtop = 50.,
    margbot = 0.,
    #zstart = 46.5, # 46.5
    #xstart = 0., # 42.
    rstart = None,
    initial = "sphere",

    # receptor params
    rec_t = 3.82e9,
    rec_p = 0.0187,
    rec_eps = 0.0,
    ka = 1e5,
    ra = 10,
)

def receptor_params(params):
    return dict(
    exclusion = False,
    walldist = 10.,
    #minsize = 0.01, # accuracy when performing reflection

    binding = True,
    eps = params.rec_eps, # margin in addition to walldist, determines re-attempting [nm]
    t = params.rec_t, # mean of exponentially distributed binding duration [ns]
    p = params.rec_p, # binding probability for one attempt
    ka = params.ka, # (bulk) association rate constant [1/Ms]
    ra = params.ra, # radius of the association zone (w/o rMolecule) [nm]
    binding_type = "proximity",

    use_force = False, # if True, t_mean = t*exp(-|F|*dx/kT)
    dx = 0.1, # width of bond energy barrier [nm]
    )
    
NAME = "rw_wei_4"
print_calculations = True

# print calculations
if print_calculations:
    phys = nanopores.Physics()
    # calculate binding probability with data from (Wei 2012)
    kon = 20.9e6 # association rate constant [1/Ms] = binding events per second
    c = 180e-9 # concentration [M = mol/l = 1000 mol/m**3]
    cmol = c * 1e3 * phys.mol # concentration [1/m**3]
    ckon = c*kon
    
    print "Average time between events (tau_on): %.2f s (from experimental data)" % (1./ckon)
    print "Number of bindings per second: %.1f (inverse of mean tau_on)" % ckon # 3.8
    
    # Smoluchowski rate equation gives number of arrivals at pore entrance per sec
    D = phys.kT / (6. * phys.pi * phys.eta * params.rMolecule * 1e-9) # [m**2/s]
    r = 6e-9 # effective radius for proteins at pore entrance [m]
    karr = 2.*phys.pi * r * D * cmol # arrival rate
    b = c * kon / karr # bindings per event
    
    print "Number of events per second: %.1f (from Smoluchowski rate equation)" % karr
    print "=> number of bindings per event: %.1f / %.1f = %.5f (= 1 - exp(-a*p) = prob of binding at least once)" % (ckon, karr, b)
    
    # solve b = 1 - exp(-ap); p = -log(1 - b)/a
    a = 0.305
    ap = -np.log(1 - b)
    p = ap/a
    print "=> a*p = -log(1 - %.5f) = %.5f" % (b, ap)
    print
    print "Average number of attempts: a = %.5f (from many simulations with dt=1, eps=0.1)" % a
    print "=> binding probability p = a*p / a = %.5f / %.5f = %.5f" % (ap, a, p)
    #receptor_params["p"] = p

def setup_rw(params):
    pore = nanopores.get_pore(**params)
    receptor = randomwalk.Ball([8.5 - 3. + 2., 0., 40.5], 0.5) # ztop 46.5    
    rw = randomwalk.RandomWalk(pore, **params)    
    rw.add_domain(receptor, **receptor_params(params))
    return rw

rw = setup_rw(params)
randomwalk.run(rw, NAME)
#rw.save(NAME)
#randomwalk.load_results(name)