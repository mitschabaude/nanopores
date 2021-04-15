# -*- coding: utf-8 -*
# (c) 2017 Gregor Mitscha-Baude
import os.path
import numpy as np
import nanopores.plots as plots
import matplotlib
matplotlib.use('Agg')
# from matplotlib import rcParams, rc
# rcParams.update({
#     "font.size" : 7,
#     "axes.titlesize" : 7,
#     #"font.family" : "sans-serif",
#     #"font.sans-serif" : ["CMU Sans Serif"],
#     "lines.linewidth" : 1,
#     "lines.markersize" : 5,
# })
import matplotlib.font_manager
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import nanopores
import nanopores.models.nanopore as nanopore_model
from nanopores.models.nanopore import Iz
import nanopores.models.randomwalk as randomwalk
from nanopores.tools import fields, statistics
fields.set_dir_mega()

colors = plots.colors


path = lambda path: os.path.join(os.path.dirname(__file__), path)

from nonlinear_least_squares import NLS
# TODO: to finish this off satisfactorily, it would be nice to infer a tau on
# histogram from exponentially distributed arrivals at the pore entrance
# according to diffusion theory

params = nanopores.user_params(
    # general params
    geoname = "wei",
    dim = 2,
    rMolecule = 3, # from 101 kDa (protein A/G/L)
    h = 5.,
    Nmax = 1e5,
    Qmol = -50., # estimate ~ protein G x2 (according to weight)
    bV = -0.2,
    dp = 26.,
    geop = dict(dp = 26.),
    posDTarget = True,

    # random walk params
    N = 20000, # number of (simultaneous) random walks
    dt = 1., # time step [ns] # .5
    walldist = 1.0, # in multiples of radius, should be >= 1
    margtop = 20.,
    margbot = 80.,
    #zstart = -46.5, # 46.5
    #xstart = 0., # 42.
    #rstart = 30.,
    initial = "bottom-disc",

    # receptor params
    tbind = 40e9, # from Lata, = 1/kd = 1/(25e-3)s [ns]
    # tbind = 286e9, from Wei, = 1/kd = 1/(3.5e-3)s [ns]
    ka = 1.5e5, # from Lata
    # ka = 3.362e8 from Wei (complicated inferred)
    # kon = 2.09e7 from Wei (simple taken)
    zreceptor = .95, # receptor location relative to pore length (1 = top)
)
##### what to do
NAME = "rw_wei_reverse_"
print_calculations = False
print_rw = False
run_test = False
compute_current = False
plot_attempt_time = False
plot_distribution = True
plot_cdf = False
voltage_dependence = False
determine_delta = False
fit_koff0 = False

##### constants
rrec = 0.5 # receptor radius
distrec = 2.75 #4. - params.rMolecule - rrec # distance of rec. center from wall
ra = distrec #params.rMolecule*(params.walldist - 1.) - rrec
dx = 0.55
kd = 25e-3
ka = 1.5e5

#### color code
color_lata = "C0" #"#0066ff"
color_wei = "#00cc00"
color_exp = "red"

def receptor_params(params):
    dx0 = params["dx"] if "dx" in params else dx
    kd0 = params["kd"] if "kd" in params else kd
    ka0 = params["ka"] if "ka" in params else ka
    tbindfromkd = 1e9/kd0
    tbind0 = params["tbind"] if "tbind" in params else tbindfromkd
    return dict(
    exclusion = False,
    walldist = 1.,
    #minsize = 0.01, # accuracy when performing reflection

    binding = True,
    t = tbind0, # mean of exponentially distributed binding duration [ns]
    #t = 1e9/kd0,
    ka = ka0, # (bulk) association rate constant [1/Ms]
    ra = ra, # radius of the association zone [nm]
    bind_type = "zone",
    collect_stats_mode = True,

    use_force = True, # if True, t_mean = t*exp(-|F|*dx/kT)
    dx = dx0, # width of bond energy barrier [nm]
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
    print('Dbulk', D)
    r = 6e-9 # effective radius for proteins at pore entrance [m] # WHY??
    r = ((params.dp - 6.)/2. - params.rMolecule)*1e-9 # should be good
    karr = 2.*phys.pi * r * D * cmol # arrival rate [1/s]
    print('karr', karr)
    b = c * kon / karr # bindings per event [1]
    
    print "Number of events per second: %.1f (from Smoluchowski rate equation)" % karr
    print "=> number of bindings per event: %.1f / %.1f = %.5f" % (ckon, karr, b)
    # (= 1 - exp(-a*p) = prob of binding at least once)
    
    # Vbind = 103 1/M
    Rbind = params.rMolecule + ra
    sin70 = np.sin(70 * np.pi / 180)
    hdist = distrec * sin70
    hseg = ra + hdist # < Rbind = ra + rMolec in our case
    Vbind = (np.pi / 3.) * hseg**2 * (3*Rbind - hseg)
    #Vbind = 4./3.*np.pi*(Rbind**3 - params.rMolecule**3) # [nm**3]
    Vbind *= (1e-8)**3 * nanopores.mol # [dm**3/mol = 1/M]
    print('Vbind', Vbind)
    cb = 0.01 # receptor concentration in binding zone [M]
    cb = 1./Vbind
    ta_ = 23.5e-9 # 1.595e-9 # avg time spent in binding zone per event [s]
    # avg no bindings = Ra * ta = (ka * cb) * ta

    ka_ = b / (cb * ta_) # implied free-solution association rate constant
    # = kon [1/Ms] * (c * V / ta) [1/s] * (1/karr) [s]
    
    print "Free solution association rate constant ka, inferred from tau on + simulations: %.4g" % ka_
    # = 3.362e08 1/Ms
    print "Simple uncorrected association rate constant in pore, from tau on: %.4g" % kon
    print "Relative correction thanks to simulation: %.1f%%" % (ka_/kon * 100,)
    # solve b = 1 - exp(-ap); p = -log(1 - b)/a
    # a = 0.305
    # ap = -np.log(1 - b)
    # p = ap/a
    # print "=> a*p = -log(1 - %.5f) = %.5f" % (b, ap)
    # print
    # print "Average number of attempts: a = %.5f (from many simulations with dt=1, eps=0.1)" % a
    # print "=> binding probability p = a*p / a = %.5f / %.5f = %.5f" % (ap, a, p)
    #receptor_params["p"] = p
    print('\n==============\n\n')

def setup_rw(params, kd=None):
    if kd:
        params["kd"] = kd
    pore = nanopores.get_pore(**params)
    rw = randomwalk.RandomWalk(pore, **params)    
    
    zrec = rw.zbot + rrec + (rw.ztop - rw.zbot - 2.*rrec)*params["zreceptor"]
    xrec = pore.radius_at(zrec) - distrec
    posrec = [xrec, 0., zrec]
    print "Receptor position: %s" % posrec
    receptor = randomwalk.Ball(posrec, rrec) # ztop 46.5
    rw.add_domain(receptor, **receptor_params(params))

    # reversed stopping
    def success(self, r, z):
        return (r**2 + (z - self.ztop)**2 > self.params.margtop**2) & (
               self.above_channel(r, z))
    def fail(self, r, z):
        return (r**2 + (z - self.zbot)**2 > self.params.margbot**2) & (
               self.below_channel(r, z))
    rw.set_stopping_criteria(success, fail)

    return rw

##### log stuff about rw
if print_rw:
    rw = randomwalk.get_rw(NAME, params, setup=setup_rw)
    print("rMolecule", params.rMolecule)
    D = rw.phys.DTargetBulk
    print("Dbulk", D)
    rporeeff = (params.dp - 6.)/2. - params.rMolecule
    print("rpore", rporeeff)
    rporeeff *= 1e-9
    # rporeeff = 6e-9 # WHY
    cmol = 180e-9 * (1e3*rw.phys.mol)
    karr = 2.*rw.phys.pi* rporeeff *D*cmol # events/s
    print("karr", karr)
    for domain in rw.domains:
        if not domain.binding or not domain.bind_type == "zone":
            continue
        print("Vbind", domain.Vbind)

    print( "percentage no attempt", 1. * np.count_nonzero(rw.attempt_times == 0) / len(rw.attempt_times) )
    print('\n==============\n\n')

##### run test rw
if run_test:
    rw = setup_rw(dict(params, bV=0., dt=10., N=500))
    randomwalk.run(rw)

##### compute current
if compute_current:
    #rw = setup_rw(params)
    #pore = nanopores.get_pore(**params)
    # Z = np.array([(rw.zbot + rw.ztop)*.5])
    # data = Iz(Z, name="current_wei_acs", calc=True, **params)
    # print(data)
    #z = (rw.zbot + rw.ztop)*0.5
    z = 41.
    print('z', z)
    x0 = None if z is None else [0., 0., z]
    params_ = dict(params, x0=x0, tol=1e-4, Nmax=2e5)
    #params_ = dict(params, x0=[0., 0., z])
    setup = nanopore_model.Setup(**params_)
    setup.geo.plot_subdomains()
    setup.geo.plot_boundaries()
    print params_
    print 'geo', setup.geo
    pb, pnps = nanopore_model.solve(setup, visualize=True)
    result = pnps.evaluate(setup.phys.CurrentPNPSDetail)
    print('result', result)
    
##### draw bindings and forces from empirical distribution
def draw_empirically(rw, N=1e8, nmax=1000, success=True, determine_ka=False):
    domain = rw.domains[1]
    N = int(N)

    Rbind = rw.params.rMolecule + domain.ra
    sin70 = np.sin(70 * np.pi / 180)
    hdist = distrec * sin70
    hseg = domain.ra + hdist # < Rbind = ra + rMolec in our case
    Vbind = (np.pi / 3.) * hseg**2 * (3*Rbind - hseg) # sphere segment volume
    Vbind *= (1e-8)**3 * nanopores.mol # [dm**3/mol = 1/M]

    ka = domain.ka
    if determine_ka:
        kon = 20.9e6 # association rate constant [1/Ms] = binding events per second
        c = 180e-9 # concentration [M = mol/l = 1000 mol/m**3]
        cmol = c * 1e3 * rw.phys.mol # concentration [1/m**3]
        print "Number of bindings per second: %.1f (inverse of mean tau_on)" % (c*kon,) # 3.8
        # Smoluchowski rate equation
        D = rw.phys.DTargetBulk
        r = ((params.dp - 6.)/2. - params.rMolecule)*1e-9
        if rw.params.initial == "bottom-disc":
            r = (rw.rbot - rw.params.rMolecule)*1e-9
        karr = 2.*rw.phys.pi * r * D * cmol # arrival rate [1/s]
        kb = c * kon / karr # bindings per event [1]
        print "Number of events per second: %.1f (from Smoluchowski rate equation)" % karr
        print "=> number of bindings per event: %.1f / %.1f = %.5f" % (c*kon, karr, kb)
        ta = 1e-9*rw.attempt_times.mean()
        ka = kb * Vbind / ta
        print "determined ka", ka

    Ra = ka / Vbind
    print "multiplying attempt times with Ra * ns = %.3g" % (Ra*1e-9,)
    # draw indices of existing random walks
    I = np.random.randint(rw.N, size=(N,))
    times = (1e-9*rw.times)[I]
    bindings = np.zeros(N, dtype=bool)
    avgbindings = (Ra*1e-9*rw.attempt_times)[I]
    bindings[avgbindings > 0] = np.random.poisson(avgbindings[avgbindings > 0])
    del avgbindings
    ibind, = np.nonzero(bindings > 0)
    n0 = len(ibind)
    n = min(n0, nmax)
    N_ = ibind[n-1] + 1
    ibind = ibind[:n]
    print "%d binding events drawn, %s used, total events used %s." % (n0, n, N_)
    print ("percentage no binding:", 1. - 1.*n / N_)
    I = I[:N_]
    times = times[:N_]
    bindings = bindings[:N_]
    
    f = np.array([f for F in rw.binding_zone_forces for f in F])
    F = np.random.choice(f, size=(n,))
    dx = 1e-9*domain.dx
    kT = rw.phys.kT
    t = domain.t * np.exp(-F*dx/kT)
    print "dwell time reduction by force:", np.mean(t)/domain.t
    bind_times = 1e-9*np.random.gamma(bindings[ibind], scale=t)
    times[ibind] += bind_times
    
    if success:
        tfail = times[rw.fail[I]]
        tsuccess = times[rw.success[I]]
        print "success", len(tsuccess), "fail", len(tfail), "percentage success", len(tsuccess)*100./(1.*len(tsuccess) + len(tfail))
        return tfail, tsuccess
    else:
        return times[ibind]

##### load tau_off histogram from source and create fake data
def tauoff_wei():
    csvfile = path("tau_off_wei.csv")
    data = np.genfromtxt(csvfile, delimiter=",")
    bins = data[:, 0]
    counts = data[:, 1]
    
    # inspection showed that there seems to be a good,
    # evenly spaced approximation to all bins except the first and last with
    # spacing 0.55, i.e. of the form (beta + 0.55*np.arange(0, N)) for some beta
    x = bins[:-1]
    N = len(x)
    # minimize norm(x - (beta + 0.55*np.arange(0, N)) w.r.t. beta
    #beta = x.mean() - 0.55*(N-1)/2.
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

###### determine tauoff from fit to exponential cdf 1 - exp(t/tauoff)
@fields.cache("wei_koff_6", default=dict(params, dx=0.55, N=10000, 
              dp=30., geop=dict(dp=30.), nmax=523, NN=1.5e8))
def fit_koff(nmax=523, NN=1.5e8, **params):
    tbind = params.pop("tbind")
    params["kd"] = 1e9/tbind
    params.pop("kd")
    dx = params.pop("dx")
    ka = params.pop("ka")
    try:
        rw = randomwalk.get_rw(NAME, params, setup=setup_rw, calc=True, kd=1e9/tbind)
    except Exception as e:
        diffs = fields.diffs(params, NAME)
        diffs = [d for d in diffs if not "bV" in d]
        print diffs
        raise e
    rw.domains[1].dx = dx
    rw.domains[1].ka = ka
    times = draw_empirically(rw, N=NN, nmax=nmax, success=False)
    bins = np.logspace(np.log10(min(times)), np.log10(max(times)), 35)
    #bins = np.logspace(-3., 2., 35)
    hist, _ = np.histogram(times, bins=bins)
    cfd = np.cumsum(hist)/float(np.sum(hist))
    t = 0.5*(bins[:-1] + bins[1:])
    tmean = times.mean()
    #toff = NLS(t, cfd, t0=tmean)
    T = statistics.Exponential(tau=None)
    T.constants['tau'] = tmean
    T.fit(times, method='cdf', log=True)
    toff = T.tau

    koff = 1./toff
    return dict(t=t, cfd=cfd, toff=toff, tmean=tmean, koff=koff)

if plot_attempt_time:
    rw = randomwalk.get_rw(NAME, params, setup=setup_rw)
    ta = rw.attempt_times
    ta = ta[ta > 0.]
    #tt = np.logspace(-.5, 2.5, 100)
    tt = np.linspace(0.25, 250., 100)
    plt.figure("attempt_times", figsize=(2.5, 1.65))
    plt.hist(ta, bins=tt, normed=True, log=True, label="Simulations")
    ta0 = ta.mean()
    #plt.plot(tt, 1./ta0 * np.exp(-tt/ta0), label="Exp. fit,\nmean %.3gns" % ta0)
    #plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel("Attempt time [ns]")
    plt.ylabel("Rel. frequency")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles[::-1], labels[::-1], frameon=False)

##### run rw in collect mode and draw bindings from empirical distributions
if plot_distribution:
    rw = randomwalk.get_rw(NAME, params, setup=setup_rw)
    ta = rw.attempt_times
    ta = ta[ta > 0.]
    #tt = np.logspace(-.5, 2.5, 100)
    #tt = np.linspace(0.25, 600., 100)
    #plt.figure("attempt_times", figsize=(4, 3))
    #plt.hist(ta, bins=tt, normed=True, log=True, label="Simulations")
    #ta0 = ta.mean()
    #plt.plot(tt, 1./ta0 * np.exp(-tt/ta0), label="Exp. fit, mean=%.3gns" % ta0)
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.xlabel("Attempt time [ns]")
    #plt.ylabel("Rel. frequency")
    #plt.legend()
    
    forces = np.array([f for F in rw.binding_zone_forces for f in F])
    plt.figure("force", figsize=(4,3))
    plt.hist(1e12*forces, bins=200, normed=True)
    plt.xlabel("Force [pN]")
    plt.ylabel("Rel. frequency")
    
    plt.figure("hist_old", figsize=(4.5,3))
    NN = 1.5e8
    fake = tauoff_wei()
    tfail, tsuccess1 = draw_empirically(rw, N=NN, nmax=len(fake))
    a, b = -6.5, 2 # log10 of plot interval
    bins = np.logspace(a, b, 40)
    _, _, gptchs = plt.hist(tsuccess1, bins=bins, color="green", log=True,
             alpha=0.6, rwidth=0.9, label=r"Translocated ($k_d$: Lata)", zorder=50)
    _, _, rptchs = plt.hist(tfail, bins=bins, color="red", log=True,
             alpha=0.6, rwidth=0.9, label=r"Did not translocate ($k_d$: Lata)")
    handler_sim1 = (gptchs[0], rptchs[0])
    
    # add histogram for kd fitted from wei
    params2 = dict(params)
    tbind = params2.pop("tbind")
    params2["kd"] = 4.5e-3 #TODO
    rw = randomwalk.get_rw(NAME, params2, setup=setup_rw)
    domain = rw.domains[1]
    #domain.ka = 1.004e+08 # 2.305e+08 #2.09e7 # hack; change ka for second histogram
    domain.initialize_binding_zone(rw)
    tfail, tsuccess2 = draw_empirically(rw, N=NN, nmax=len(fake), determine_ka=True)
    _, _, gptchs = plt.hist(tsuccess2, bins=bins, color="green", log=True,
                            histtype="step", linestyle="--", alpha=0.6,
                            rwidth=0.9, label=r"Translocated ($k_d$: Wei)", zorder=200)
    _, _, rptchs = plt.hist(tfail, bins=bins, color="red", log=True,
                            histtype="step", linestyle="--", alpha=0.6,
                            rwidth=0.9, label=r"Did not translocate ($k_d$: Wei)")
    handler_sim2 = (gptchs[0], rptchs[0])
    
    _, _, ptchs = plt.hist(fake, bins=bins, histtype="step", log=True,
             color="orange", label="Experiments", zorder=100)
    handler_exp = ptchs[0]
    
    plt.xscale("log")
    #plt.yscale("log")
    plt.ylabel("Count")
    plt.xlabel(r"$\tau$ off [s]")
    plt.ylim(ymin=1.)
    #plt.xlim(xmax=1e4)
    #plt.legend()
    # plt.legend([handler_exp, handler_sim1, handler_sim2],
    #            ["Experiments", r"Sim. ($k_d$ from Lata)", r"Sim. ($k_d$ from Wei)"],
    # #plt.legend([handler_exp, handler_sim1],
    # #           ["Experiments (Wei et al.)", r"Sim. ($k_d$ from Lata)"],
    #             handler_map={tuple: HandlerTuple(ndivide=None)},
    #             frameon=False)
    #scatterpoints=1, numpoints=1, 
    
    # simpler hist with better color code and visibility
    tsuccess_lata = tsuccess1#[tsuccess1 > 1e-4]
    tsuccess_wei = tsuccess2#[tsuccess2 > 1e-4]
    plt.figure("hist_all", figsize=(2.75, 1.83333333333))
    plt.hist(tsuccess_lata, bins=bins, color=color_lata, log=True,
             alpha=0.8, rwidth=0.9, label=r"Sim. ($k_a$, $k_d$ from Lata)", zorder=50)
    
    plt.hist(tsuccess_wei, bins=bins, color=color_wei, log=True,
             #histtype="step", linestyle="--", 
             alpha=0.5, rwidth=0.9, label=r"Sim. ($k_a$, $k_d$ from Wei)", zorder=90)
    plt.hist(fake, bins=bins, histtype="step", log=True,
             linewidth=1.75,
             color=color_exp, label="Experiment", zorder=100)
    plt.xscale("log")
    #plt.yscale("log")
    plt.ylabel("Count")
    plt.xlabel(r"$\tau_\mathrm{off}$ [ms]")
    plt.ylim(ymin=1.)
    ax = plt.gca()
    ax.set_xticks([1e-6, 1e-3, 1.])
    ax.set_xticks([1e-5, 1e-4, 1e-2, 1e-2, 1e-1, 1e1, 1e2], minor=True)
    ax.set_xticklabels(["$\mathregular{10^{-3}}$", "1", "$\mathregular{10^3}$"])
    ax.set_xticklabels([], minor=True)
    #plt.xlim(xmin=.3e-6, xmax=1e2)
    #plt.xlim(xmin=0.2e-4, xmax=0.9e2)
    plt.legend(loc="best", frameon=False)

###### reproduce cumulative tauoff plot with fits and different bV
voltages = [-0.2, -0.25, -0.3, -0.35][::-1]
colorsList = ["k", "r", "b", "g"][::-1]
zrecs = [.90, .95, .99]
N = 10000
newparams = dict(N=N, dp=30., geop=dict(dp=30.), margbot=70.)

if plot_cdf:
    plt.figure("bV_tauoff", figsize=(2.5, 1.9))
    for i, v in enumerate(voltages):
        data = fit_koff(bV=v, zreceptor=.95, dx=dx, **newparams)
        tt = np.logspace(np.log10(np.min(data.t)*.5), np.log10(np.max(data.t)*2), 100)
        
        lines = plt.semilogx(tt, 1. - np.exp(-tt/data.toff), color=colorsList[i])
        plt.semilogx(data.t, data.cfd, "v", color=lines[0].get_color(),
                     label="%d mV" % (1000*abs(v)))
        print "koff", data.koff
    plt.xlim(1e-4, 1e1)
    plt.xlabel(r"$\tau$ off [s]")
    plt.ylabel("Cumulative probability")
    plt.legend(frameon=False)
    
###### regression to quantify bV-koff relationship
def regression(bV, koff):
    "find coefficients in relationship koff = koff0 * exp(a*bV)"
    X = np.column_stack([bV, np.ones(len(bV))])
    y = np.log(koff)
    a, b = tuple(np.dot(np.linalg.inv(np.dot(X.T, X)), np.dot(X.T, y)))
    return a, np.exp(b)

######
def koff0(kd, **params):
    return fit_koff(name="wei_koff_6", bV=0., tbind=1e9/kd, margbot=70., **params).koff

if fit_koff0:
    plt.figure("koff0", figsize=(2.5, 1.9))
    kd = np.array([3., 4, 4.25, 4.5, 4.75, 5, 6.])
    ko = np.array([1e3*koff0(k*1e-3, ka=1e8, dx=0.55, NN=1e8, nmax=10000) for k in kd])
    c = ko.mean()/kd.mean()
    print "F0 = %.3f pN" % (1e12*np.log(c)/(1e-9*0.55)*nanopores.kT)
    # F0 = 0.184 pN
    plt.axhline(y=4.5, linestyle="-", color="C0", label="Wei et al.")
    plt.plot(kd, ko, "oC1", label="Simulations")
    #plt.plot(kd, c*kd, ":C1", label=r"Fit to $k_{off}^{V=0}$ = C*$k_d$")
    plt.plot(kd, kd, ":C1", label=r"$k_{off}^{V=0}$ = $k_d$")
    plt.ylim(ymin=2.5, ymax=6.5)
    #plt.xlim(2.9, 4.6)
    #plt.fill_between(plt.xlim(), [4.5 - 0.6]*2, [4.5 + 0.6]*2, alpha=0.5, color="C1")
    
    #plt.xticks(kd)
    #plt.yticks(kd)
    #plt.xlim(plt.ylim())
    #plt.axis("equal")
    plt.xlabel(r"Bulk dissociation constant $k_d$ [10$^{-3}$/(Ms)]")
    plt.ylabel(r"$k_{off}^{V=0}$ [10$^{-3}$/(Ms)]")
    plt.legend(frameon=False, loc="upper left")

###### recreate voltage-dependent plot of koff
if voltage_dependence:
    # get experimental data
    plt.figure("koff", figsize=(2.2, 1.83333333))
    data = np.genfromtxt(path("koff.csv"), delimiter=",")
    v = data[:, 0]
    koff = data[:, 1]
    c0, k0 = regression(np.abs(v), koff)
    vv = np.linspace(0., 370., 10)
    plt.plot(vv, k0 * np.exp(c0*vv), "-r", lw=1.75)
    
    v = np.array([-0., -0.05, -0.1, -0.15, -0.2, -0.25, -0.3, -0.35])
    #v = np.array([-0.2, -0.25, -0.3, -0.35])
    mv = np.abs(v)*1e3
    z = 0.95
    #dx = 0.55
    koff = [fit_koff(bV=V, zreceptor=z, dx=dx, **newparams).koff for V in v]
    c1, k1 = regression(mv[-4:], koff[-4:])
    plt.plot(mv, koff, "v", markersize=7, label=r"Sim. ($k_d$ from Lata)",
             color=color_lata)
    plt.plot(vv, k1 * np.exp(c1*vv), "-", color=color_lata)
             
    # and again with different kd
    kd = 4.5e-3
    koff1 = [fit_koff(bV=V, zreceptor=z, dx=dx,
                     tbind=1e9/kd, **newparams).koff for V in v]
    c2, k2 = regression(mv[-4:], koff1[-4:])
    plt.plot(mv, koff1, "o", markersize=7, label=r"Sim. ($k_d$ from Wei)",
             color=color_wei)
             #mfc="None", mec=color_wei)
    #plt.plot(vv, k2 * np.exp(c2*vv), ":", color="#990000")
    
    v = data[:, 0]
    koff2 = data[:, 1]
    plt.plot(v, koff2, "s", mfc="None", mec=color_exp,
             markersize=6, mew=1.75, label="Experiment")

    plt.yscale("log")
    plt.ylim(ymax=.9e3)
    plt.xlabel("Voltage [mV]")
    plt.ylabel(r"$k_\mathrm{off}$ [1/s]")
    plt.legend(frameon=False, loc="upper left")
    
    plt.figure("koff_simple", figsize=(1.7, 1.6))
    plt.plot(mv, koff1, "o", markersize=7, label=r"Simulation", color=color_wei)
    plt.plot(v, koff2, "s", markersize=6, mew=1.75,
             label=r"Experiment", mec=color_exp, mfc="None")
    plt.yscale("log")
    #plt.xlabel("Voltage [mV]")
    #plt.ylabel("k off [1/s]")
    plt.ylabel(ur"$\log(k_\mathrm{off})$")
    plt.xlabel("Voltage")

    plt.tick_params(
        axis="both",          # changes apply to the x-axis
        which="both",      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        left=False,
        right=False,
        labelleft=False,
        labelbottom=False)
    plt.legend(frameon=False)
    plt.title("Reaction kinetics")
    
    
###### read koff-bV dependence from wei data
koff0 = np.array([])
coeff = np.array([])
for i in range(1, 6):
    data = np.genfromtxt(path("koff%d.csv" %i), delimiter=",")
    voltages = data[:, 0]*1e-3
    koff = data[:, 1]
    
    c, k = regression(np.abs(voltages), koff)
    coeff = np.append(coeff, c)
    koff0 = np.append(koff0, k)

cdxall_exp = coeff
cdx_exp = coeff.mean()
vdx_exp = coeff.std()

###### plt determination of bond rupture length from wei data and simulations
if determine_delta:
    voltages = [-0.2, -0.25, -0.3, -0.35]
    zrecs = [.91, .93, .95, .97]
    dxtest = 0.4
    dx = dxtest
    koff0 = np.array([])
    coeff = np.array([])
    for z in zrecs:
        for v, koff in nanopores.collect(voltages):
            data = fit_koff(bV=v, zreceptor=z, dx=dx, NN=4e8, nmax=2000, **dict(newparams,  margbot=80.))
            koff.new = data.koff
        c, k = regression(np.abs(voltages), koff)
        coeff = np.append(coeff, c)
        koff0 = np.append(koff0, k)
    
    cdxtest_sim = coeff.mean()
    
    dx0 = cdx_exp / (cdxtest_sim / dxtest)
    print "inferred dx:", dx0
    
    dxs = [0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8] #2., 3., 4., 4.5, 4.79, 5.1, 5.5, 6., 7.]
    cdx = []
    cdxstd = []
    cdxall = []
    for dx in dxs:
        #print "dx", dx    
        koff0 = np.array([])
        coeff = np.array([])
        for z in zrecs:
            for v, koff in nanopores.collect(voltages):
                data = fit_koff(NN=1e8, ka=1e8, nmax=10000, bV=v, zreceptor=z, dx=dx, **dict(newparams,  margbot=80.))
                koff.new = data.koff
            c, k = regression(np.abs(voltages), koff)
            coeff = np.append(coeff, c)
            koff0 = np.append(koff0, k)
        cdx.append(coeff.mean())
        cdxall.append(coeff)
        cdxstd.append(coeff.std())
        #print "c*dx %.3g +- %.3g" % (coeff.mean(), coeff.std())
        #print "koff0 %.3g +- %.3g" % (koff0.mean(), koff0.std())
    def fplot(cdx, dx):
        return cdx
        
    dx = np.array(dxs)
    cdx = np.array(cdx)
    cdxall = np.array(cdxall)
    cdxstd = np.array(cdxstd)
    
    dxx = np.linspace(dx[0], dx[-1], 100)
    cdx_exp = np.ones(len(dxx))*cdx_exp
    vdx_exp = np.ones(len(dxx))*vdx_exp
    
    plt.figure("delta", figsize=(2.5, 1.9))
    plt.plot(dxx, fplot(cdx_exp, dxx), "-", label="Wei et al.")
    for i in range(5):
        plt.plot(dxx, np.ones(len(dxx))*cdxall_exp[i], "-", color="C0", alpha=0.5)
    #plt.plot(dxx, fplot(cdx_exp - vdx_exp, dxx), "-", color="C1")
    #plt.plot(dxx, fplot(cdx_exp + vdx_exp, dxx), "-", color="C1")
    #plt.fill_between(dxx, fplot(cdx_exp - vdx_exp, dxx),
    #                     fplot(cdx_exp + vdx_exp, dxx), alpha=0.5)
    
    #plt.plot(dx, fplot(cdx, dx), "o", label="Simulation")
    plt.plot(dx, fplot(cdx, dx), "o", label="Simulations", color="C1")
    for i in range(cdxall.shape[1]):
        plt.plot(dx, fplot(cdxall[:, i], dx), "o", color="C1", alpha=0.5)
    #plt.fill_between(dx, fplot(cdx - cdxstd, dx), fplot(cdx + cdxstd, dx), alpha=0.5)
    
    plt.annotate(r"$\delta$=0.55nm", (0.55, cdxall[3, 0] - .6),
                 xytext=(0.55 - .03, cdxall[3, 0] - 5.), color="C1",
                 arrowprops=dict(arrowstyle="->", color="C1"))
    
    plt.xlabel(r"Bond rupture length $\delta$ [nm]")
    plt.ylabel(r"$\alpha$ [1/V]")
    plt.legend(loc="upper left", frameon=False)

import folders
nanopores.savefigs("tau_off4", folders.FIGDIR + "/wei", ending=".pdf")
#nanopores.savefigs("tau_off4", folders.FIGDIR_HOWORKA + "/wei", ending=".pdf")
#nanopores.savefigs("tau_off2", folders.FIGDIR + "/wei", ending=".eps")
