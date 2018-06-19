# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import nanopores
import nanopores.models.randomwalk as randomwalk
from nanopores.tools import fields
fields.set_dir_mega()

from nonlinear_least_squares import NLS
# TODO: to finish this off satisfactorily, it would be nice to infer a tau on
# histogram from exponentially distributed arrivals at the pore entrance
# according to diffusion theory

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
    tbind = 40e9, # from Lata, = 1/kd = 1/(25e-3)s [ns]
    # tbind = 286e9, from Wei, = 1/kd = 1/(3.5e-3)s [ns]
    ka = 1.5e5,
    zreceptor = .95, # receptor location relative to pore length (1 = top)
)
##### what to do
NAME = "rw_wei_"
print_calculations = False
run_test = False
plot_distribution = False
plot_cdf = False
voltage_dependence = True
determine_delta = False
fit_koff0 = False

##### constants
rrec = 0.5 # receptor radius
distrec = 4. - params.rMolecule - rrec # distance of rec. center from wall
ra = distrec #params.rMolecule*(params.walldist - 1.) - rrec
dx = 5.5
kd = 25e-3

#### color code
color_lata = "C0" #"#0066ff"
color_wei = "#00cc00"
color_exp = "red"

def receptor_params(params):
    dx0 = params["dx"] if "dx" in params else dx
    kd0 = params["kd"] if "kd" in params else kd
    tbindfromkd = 1e9/kd0
    tbind0 = params["tbind"] if "tbind" in params else tbindfromkd
    return dict(
    exclusion = False,
    walldist = 1.,
    #minsize = 0.01, # accuracy when performing reflection

    binding = True,
    t = tbind0, # mean of exponentially distributed binding duration [ns]
    #t = 1e9/kd0,
    ka = params["ka"], # (bulk) association rate constant [1/Ms]
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
def draw_empirically(rw, N=1e8, nmax=1000, success=True):
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
    
    if success:
        tfail = times[rw.fail[I]]
        tsuccess = times[rw.success[I]]
        return tfail, tsuccess
    else:
        return times[ibind]

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
@fields.cache("wei_koff_2", default=dict(params, dx=5.5, N=10000, 
              dp=30., geop=dict(dp=30.), nmax=523, NN=4e8))
def fit_koff(nmax=523, NN=4e8, **params):
    tbind = params.pop("tbind")
    params["kd"] = 1e9/tbind
    dx = params.pop("dx")
    rw = randomwalk.get_rw(NAME, params, setup=setup_rw, calc=True)
    rw.domains[1].dx = dx
    times = draw_empirically(rw, N=NN, nmax=nmax, success=False)
    bins = np.logspace(np.log10(min(times)), np.log10(max(times)), 35)
    #bins = np.logspace(-3., 2., 35)
    hist, _ = np.histogram(times, bins=bins)
    cfd = np.cumsum(hist)/float(np.sum(hist))
    t = 0.5*(bins[:-1] + bins[1:])
    tmean = times.mean()
    toff = NLS(t, cfd, t0=tmean)
    koff = 1./toff
    return dict(t=t, cfd=cfd, toff=toff, tmean=tmean, koff=koff)

##### run rw in collect mode and draw bindings from empirical distributions
if plot_distribution:
    rw = randomwalk.get_rw(NAME, params, setup=setup_rw)
    ta = rw.attempt_times
    ta = ta[ta > 0.]
    #tt = np.logspace(-.5, 2.5, 100)
    tt = np.linspace(0.25, 200., 100)
    plt.figure("attempt_times", figsize=(4,3))
    plt.hist(ta, bins=tt, normed=True, log=True, label="Simulations")
    ta0 = ta.mean()
    plt.plot(tt, 1./ta0 * np.exp(-tt/ta0), label="Exp. fit, mean=%.3gns" % ta0)
    #plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel("Attempt time [ns]")
    plt.ylabel("Rel. frequency")
    plt.legend()
    
    forces = np.array([f for F in rw.binding_zone_forces for f in F])
    plt.figure("force", figsize=(4,3))
    plt.hist(1e12*forces, bins=200, normed=True)
    plt.xlabel("Force [pN]")
    plt.ylabel("Rel. frequency")
    
    plt.figure("hist_old", figsize=(4.5,3))
    NN = 3e8
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
    params2["kd"] = 3.5e-3
    rw = randomwalk.get_rw(NAME, params2, setup=setup_rw)
    tfail, tsuccess2 = draw_empirically(rw, N=NN, nmax=len(fake))
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
    plt.legend([handler_exp, handler_sim1, handler_sim2],
               ["Experiments", r"Sim. ($k_d$ from Lata)", r"Sim. ($k_d$ from Wei)"],
    #plt.legend([handler_exp, handler_sim1],
    #           ["Experiments (Wei et al.)", r"Sim. ($k_d$ from Lata)"],
                handler_map={tuple: HandlerTuple(ndivide=None)},
                frameon=False)
    #scatterpoints=1, numpoints=1, 
    
    # simpler hist with better color code and visibility
    tsuccess_lata = tsuccess1#[tsuccess1 > 1e-4]
    tsuccess_wei = tsuccess2#[tsuccess2 > 1e-4]
    plt.figure("hist_all", figsize=(4.5,3))
    plt.hist(tsuccess_lata, bins=bins, color=color_lata, log=True,
             alpha=0.8, rwidth=0.9, label=r"Sim. ($k_d$ from Lata)", zorder=50)
    
    plt.hist(tsuccess_wei, bins=bins, color=color_wei, log=True,
             #histtype="step", linestyle="--", 
             alpha=0.5, rwidth=0.9, label=r"Sim. ($k_d$ from Wei)", zorder=90)
    plt.hist(fake, bins=bins, histtype="step", log=True,
             linewidth=3,
             color=color_exp, label="Experiment", zorder=100)
    plt.xscale("log")
    #plt.yscale("log")
    plt.ylabel("Count")
    plt.xlabel(r"$\tau$ off [s]")
    plt.ylim(ymin=1.)
    #plt.xlim(xmin=.3e-6, xmax=1e2)
    #plt.xlim(xmin=0.2e-4, xmax=0.9e2)
    plt.legend(loc="best", frameon=False)

###### reproduce cumulative tauoff plot with fits and different bV
voltages = [-0.2, -0.25, -0.3, -0.35][::-1]
colors = ["k", "r", "b", "g"][::-1]
zrecs = [.90, .95, .99]
N = 10000
newparams = dict(N=N, dp=30., geop=dict(dp=30.))

if plot_cdf:
    plt.figure("bV_tauoff", figsize=(4, 3))
    for i, v in enumerate(voltages):
        data = fit_koff(bV=v, zreceptor=.95, dx=dx, **newparams)
        tt = np.logspace(-3., 2., 100)
        
        lines = plt.semilogx(tt, 1. - np.exp(-tt/data.toff), color=colors[i])
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
    return fit_koff(name="wei_koff_3", bV=0., tbind=1e9/kd, **params).koff

if fit_koff0:
    plt.figure("koff0", figsize=(4, 3))
    kd = np.array([2., 3, 3.25, 3.5, 3.75, 4, 5.])
    ko = np.array([1e3*koff0(k*1e-3, NN=5e8, nmax=10000) for k in kd])
    c = ko.mean()/kd.mean()
    print "F0 = %.3f pN" % (1e12*np.log(c)/(1e-9*5.5)*nanopores.kT)
    # F0 = 0.184 pN
    plt.axhline(y=4.5, linestyle="-", color="C0", label="Wei et al.")
    plt.plot(kd, ko, "oC1", label="Simulations")
    plt.plot(kd, c*kd, ":C1", label=r"Fit to $k_{off}^{V=0}$ = C*$k_d$")
    plt.plot(kd, kd, ":C2", label=r"$k_{off}^{V=0}$ = $k_d$")
    plt.ylim(ymin=2.2, ymax=7.5)
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
    plt.figure("koff", figsize=(3.6, 3))
    data = np.genfromtxt("koff.csv", delimiter=",")
    v = data[:, 0]
    koff = data[:, 1]
    c0, k0 = regression(np.abs(v), koff)
    vv = np.linspace(0., 400., 10)
    plt.plot(vv, k0 * np.exp(c0*vv), "-r", lw=2)
    
    v = np.array([-0., -0.05, -0.1, -0.15, -0.2, -0.25, -0.3, -0.35])
    mv = np.abs(v)*1e3
    z = 0.95
    dx = 5.6
    koff = [fit_koff(bV=V, zreceptor=z, dx=dx, **newparams).koff for V in v]
    c1, k1 = regression(mv[-4:], koff[-4:])
    plt.plot(mv, koff, "v", markersize=10, label=r"Sim. ($k_d$ from Lata)",
             color=color_lata)
    plt.plot(vv, k1 * np.exp(c1*vv), "-", color=color_lata)
             
    # and again with different kd
    kd = 3.5e-3
    koff1 = [fit_koff(bV=V, zreceptor=z, dx=dx,
                     tbind=1e9/kd, **newparams).koff for V in v]
    c2, k2 = regression(mv[-4:], koff1[-4:])
    plt.plot(mv, koff1, "o", markersize=10, label=r"Sim. ($k_d$ from Wei)",
             color=color_wei)
             #mfc="None", mec=color_wei)
    #plt.plot(vv, k2 * np.exp(c2*vv), ":", color="#990000")
    
    v = data[:, 0]
    koff2 = data[:, 1]
    plt.plot(v, koff2, "s", mfc="None", mec=color_exp,
             markersize=8, mew=3, label="Experiment")

    plt.yscale("log")
    plt.ylim(ymax=.9e3)
    plt.xlabel("Voltage [mV]")
    plt.ylabel("k off [1/s]")
    plt.legend(frameon=False, loc="upper left")
    
    plt.figure("koff_simple", figsize=(2.7, 2.3))
    plt.plot(mv, koff, "v", markersize=10, label=r"Simulation", color="C0")
    plt.plot(v, koff2, "s", markersize=8, label=r"Experiment", color="red")
    plt.xlabel("Voltage [mV]")
    plt.ylabel("k off [1/s]")
    plt.yscale("log")
#    plt.tick_params(
#        axis="both",          # changes apply to the x-axis
#        which="both",      # both major and minor ticks are affected
#        bottom=False,      # ticks along the bottom edge are off
#        top=False,         # ticks along the top edge are off
#        left=False,
#        right=False,
#        labelleft = False,
#        labelbottom = False
#        ) # labels along the bottom edge are off
    plt.legend(frameon=False)
    
    
###### read koff-bV dependence from wei data
koff0 = np.array([])
coeff = np.array([])
for i in range(1, 6):
    data = np.genfromtxt("koff%d.csv" %i, delimiter=",")
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
    zrecs = [.90, .95, .99]
    dxtest = 5.
    dx = dxtest
    koff0 = np.array([])
    coeff = np.array([])
    for z in zrecs:
        for v, koff in nanopores.collect(voltages):
            data = fit_koff(bV=v, zreceptor=z, dx=dx, **newparams)
            koff.new = data.koff
        c, k = regression(np.abs(voltages), koff)
        coeff = np.append(coeff, c)
        koff0 = np.append(koff0, k)
    
    cdxtest_sim = coeff.mean()
    
    dx0 = cdx_exp / (cdxtest_sim / dxtest)
    print "inferred dx:", dx0
    
    dxs = [2., 3., 4., 5., 5.5, 6., 7., 8., 9.]
    cdx = []
    cdxstd = []
    cdxall = []
    for dx in dxs:
        #print "dx", dx    
        koff0 = np.array([])
        coeff = np.array([])
        for z in zrecs:
            for v, koff in nanopores.collect(voltages):
                data = fit_koff(bV=v, zreceptor=z, dx=dx, **newparams)
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
    
    plt.figure("delta", figsize=(4, 3))
    plt.plot(dxx, fplot(cdx_exp, dxx), "-", label="Wei et al.")
    for i in range(5):
        plt.plot(dxx, np.ones(len(dxx))*cdxall_exp[i], "-", color="C0", alpha=0.5)
    #plt.plot(dxx, fplot(cdx_exp - vdx_exp, dxx), "-", color="C1")
    #plt.plot(dxx, fplot(cdx_exp + vdx_exp, dxx), "-", color="C1")
    #plt.fill_between(dxx, fplot(cdx_exp - vdx_exp, dxx),
    #                     fplot(cdx_exp + vdx_exp, dxx), alpha=0.5)
    
    #plt.plot(dx, fplot(cdx, dx), "o", label="Simulation")
    plt.plot(dx, fplot(cdx, dx), "o", label="Simulations", color="C1")
    for i in (0, 1, 2):
        plt.plot(dx, fplot(cdxall[:, i], dx), "o", color="C1", alpha=0.5)
    #plt.fill_between(dx, fplot(cdx - cdxstd, dx), fplot(cdx + cdxstd, dx), alpha=0.5)
    
    plt.annotate(r"$\delta$=5.5nm", (5.5, cdxall[4, 0] - 1.),
                 xytext=(5.5 - .79, cdxall[4, 0] - 8.), color="C1",
                 arrowprops=dict(arrowstyle="->", color="C1"))
    
    plt.xlabel(r"Bond rupture length $\delta$ [nm]")
    plt.ylabel(r"$\alpha$ [1/V]")
    plt.legend(loc="upper left", frameon=False)

import folders
nanopores.savefigs("tau_off2", folders.FIGDIR + "/wei", ending=".pdf")
