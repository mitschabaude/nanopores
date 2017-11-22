# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as  plt
import nanopores
import nanopores.models.randomwalk as randomwalk
from nanopores.tools import fields
fields.set_dir_mega()
# TODO: fit bond rupture length to wei data

params = nanopores.user_params(
    # general params
    geoname = "wei",
    dim = 2,
    rMolecule = 1.25, # 6.
    h = 5.,
    Nmax = 1e5,
    Qmol = 2., #15.,
    bV = -0.2,
    dp = 26.,
    geop = dict(dp = 26.),
    posDTarget = True,

    # random walk params
    N = 100000, # number of (simultaneous) random walks
    dt = .5, # time step [ns]
    walldist = 2., # in multiples of radius, should be >= 1
    margtop = 60.,
    margbot = 0.,
    #zstart = 46.5, # 46.5
    #xstart = 0., # 42.
    rstart = 30,
    initial = "sphere",

    # receptor params
    tbind = 40e9, # = 1/kd = 1/(25e-3)s [ns]
    ka = 1.5e5,
    zreceptor = .95, # receptor location relative to pore length (1 = top)
)
##### what to do
NAME = "rw_wei_"
print_calculations = False
run_test = False
plot_distribution = False

##### constants
rrec = 0.5 # receptor radius
distrec = 4. - params.rMolecule - rrec # distance of rec. center from wall
ra = distrec #params.rMolecule*(params.walldist - 1.) - rrec

def receptor_params(params):
    return dict(
    exclusion = False,
    walldist = 1.,
    #minsize = 0.01, # accuracy when performing reflection

    binding = True,
    t = params.tbind, # mean of exponentially distributed binding duration [ns]
    ka = params.ka, # (bulk) association rate constant [1/Ms]
    ra = ra, # radius of the association zone [nm]
    bind_type = "zone",
    collect_stats_mode = True,

    use_force = True, # if True, t_mean = t*exp(-|F|*dx/kT)
    dx = 3., # width of bond energy barrier [nm]
    )

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
    rw = randomwalk.RandomWalk(pore, **params)    
    
    zrec = rw.zbot + rrec + (rw.ztop - rw.zbot - 2.*rrec)*params["zreceptor"]
    xrec = pore.radius_at(zrec) - distrec
    posrec = [xrec, 0., zrec]
    print "Receptor position: %s" % posrec
    receptor = randomwalk.Ball(posrec, rrec) # ztop 46.5
    rw.add_domain(receptor, **receptor_params(params))
    return rw

##### run test rw
if run_test:
    rw = setup_rw(params)
    randomwalk.run(rw)
    
##### draw bindings and forces from empirical distribution
def draw_empirically(rw, N=1e8, nmax=1000):
    self = rw.domains[1]
    N = int(N)
    ka = self.kbind
    # draw indices of existing random walks
    I = np.random.randint(rw.N, size=(N,))
    times = (1e-9*rw.times)[I]
    bindings = np.zeros(N, dtype=bool)
    avgbindings = (ka*rw.attempt_times)[I]
    bindings[avgbindings > 0] = np.random.poisson(avgbindings[avgbindings > 0])
    del avgbindings
    ibind, = np.nonzero(bindings > 0)
    n0 = len(ibind)
    n = min(n0, nmax)
    ibind = ibind[:n]
    print "%d binding events drawn, %s used." % (n0, n)
    
    f = np.array([f for F in rw.binding_zone_forces for f in F])
    F = np.random.choice(f, size=(n,))
    dx = 1e-9*self.dx
    kT = rw.phys.kT
    t = self.t * np.exp(-F*dx/kT)
    print "dwell time reduction by force:", np.mean(t)/self.t
    bind_times = 1e-9*np.random.gamma(bindings[ibind], scale=t)
    times[ibind] += bind_times
    
    tfail = times[rw.fail[I]]
    tsuccess = times[rw.success[I]]
    return tfail, tsuccess

##### load tau_off histogram from source and create fake data
def tauoff_wei():
    csvfile = "tau_off_wei.csv"
    data = np.genfromtxt(csvfile, delimiter=",")
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
    
    # TODO: need better experimental data => webtool
    # create fake data samples that reproduce the histogram
    fake = np.array([])
    
    frac = 1.
    while int(counts[0]*frac) > 1:
        frac /= 2.
        a, b = bins[1]*frac, bins[1]*2*frac
        sample = a*(b/a)**(np.random.rand(int(counts[0]*frac)))
        fake = np.append(fake, sample)
        #print "frac", frac
    
    for i in range(1, N):
        a, b = bins[i], bins[i+1]
        sample = a*(b/a)**(np.random.rand(counts[i]))
        fake = np.append(fake, sample)
        
    print len(fake), "events loaded from experimental data."
    return fake

##### run rw in collect mode and draw bindings from empirical distributions
if plot_distribution:
    rw = randomwalk.get_rw(NAME, params, setup=setup_rw)
    ta = rw.attempt_times
    ta = ta[ta > 0.]
    #tt = np.logspace(-.5, 2.5, 100)
    tt = np.linspace(0.25, 200., 100)
    plt.figure("attempt_times")
    plt.hist(ta, bins=tt, normed=True, label="Simulations")
    ta0 = ta.mean()
    plt.plot(tt, 1./ta0 * np.exp(-tt/ta0), label="Exp. fit, mean=%.3gns" % ta0)
    #plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Attempt time [ns]")
    plt.ylabel("Rel. frequency")
    plt.legend()
    
    forces = np.array([f for F in rw.binding_zone_forces for f in F])
    plt.figure("force")
    plt.hist(1e12*forces, bins=200, normed=True)
    plt.xlabel("Force [pN]")
    plt.ylabel("Rel. frequency")
    
    plt.figure("hist")
    fake = tauoff_wei()
    tfail, tsuccess = draw_empirically(rw, N=3e8, nmax=len(fake))
    a, b = -6.5, 3 # log10 of plot interval
    bins = np.logspace(a, b, 40)
    plt.hist(tsuccess, bins=bins, color="green", alpha=0.6, rwidth=0.9, label="Translocated", zorder=50)
    plt.hist(tfail, bins=bins, color="red", alpha=0.6, rwidth=0.9, label="Did not translocate")
    plt.hist(fake, bins=bins, histtype="step", color="orange", label="Wei et al.", zorder=100)
    plt.xscale("log")
    plt.yscale("log")
    plt.ylabel("Count")
    plt.xlabel(r"$\tau$ off [s]")
    plt.ylim(ymin=1.)
    plt.legend()
    
# determine tauoff from fit to exponential cdf 1 - exp(t/tauoff)
voltages = [-0.2, -0.25, -0.3, -0.35]
zrecs = [.90, .95, .99]
N = 10000
params.update(N=N, dp=30., geop=dict(dp=30.))
for v in voltages:
    for z in zrecs:
        params.update(bV=v, zreceptor=z)
        rw = randomwalk.get_rw(NAME, params, setup=setup_rw)
    
# recreate voltage-dependent plot of tauoff
#params.update(
#        )
  
import folders
nanopores.savefigs("tau_off2", folders.FIGDIR + "/wei", (4, 3), ending=".pdf")
