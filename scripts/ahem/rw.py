# -*- coding: utf-8 -*-
# (c) 2017 Gregor Mitscha-Baude
"script to generate all figures for exit-time paper"
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm

import nanopores
from nanopores.models.randomwalk import (get_pore, RandomWalk, run,
                                         load_results, get_rw, setup_default)
#from nanopores.models.randomwalk import get_results
from nanopores.models.nanopore import Iz, get_active_params
from nanopores.tools import fields
from nanopores import Params
#fields.set_dir_dropbox()
fields.set_dir_mega()
#fields.set_dir(fields.HOME + "/code/nanopores/fields")
FIGDIR = nanopores.dirnames.DROPBOX_FIGDIR

########### PARAMETERS ###########
# TODO: confirm rMolecule and Qmol !!!
# goal: dt=0.01, T=10000, which means up to 1e6 steps
# take ~3.6s per 1000 samples, or ~10h for 100000 samples
params = nanopores.Params(
    # simulation params
    geoname = "alphahem",
    dim = 2,
    rMolecule = .5,
    h = 1.,
    Nmax = 1e5,
    Qmol = -2.,
    bV = 0.5,
    posDTarget = False,
    R=21, Hbot=21, Htop=21,
    geop = dict(R=21, Hbot=21, Htop=21),
    ahemqs = None,
    ahemuniformqs = False,
    
    # random walk params
    N = 10000, # number of (simultaneous) random walks
    dt = 0.01, # time step [ns]
    walldist = 1., # in multiples of radius, should be >= 1
    rstart = 1.,
    zstart = 1.,
    initial = "disc",  # oder "sphere"
    
    # binding params
    t_bind = 1000.,
    p_bind = 1.,
    eps_bind = 0.1,
    
    # stopping criterion: max time (w/o binding) and radius
    Tmax = 10000.,
    Rmax = 20.,
    zstop = -1.28,
)
NAME = "rw_exittime"

########### WHAT TO DO  ###########
plot_streamlines = False
run_test = False
plot_rw_results = True
do_calculations = False
create_current_trace = False

########### SET UP RANDOM WALK  ###########
def setup_rw(params):
    pore = get_pore(**params)
    rw = RandomWalk(pore, **params)
    rw.add_wall_binding(t=params.t_bind, p=params.p_bind, eps=params.eps_bind)
    
    # define non-standard stopping criteria
    Tmax = params.Tmax
    Rmax = params.Rmax
    
    def success(self, r, z):
        return self.in_channel(r, z) & (z <= params.zstop)
    
    def fail(self, r, z):
        if self.t > Tmax:
            return np.full(r.shape, True, dtype=bool)
        toolong = (self.times[self.alive] + self.bind_times[self.alive]) > 5e6
        toofar = r**2 + z**2 > Rmax**2
        return toolong | toofar
    
    rw.set_stopping_criteria(success, fail)
    return rw

########### STREAMLINE PLOT  ###########
if plot_streamlines:
    rw = setup_rw(params)
    rw.plot_streamlines(both=True, R=10, Hbot=15, Htop=10,
                        maxvalue=1e-10, figsize=(5, 5))
    plt.figure()

########### RUN A TEST RANDOM WALK ###########
if run_test:
    rw = setup_default(Params(params, N=100))
    #setup_rw()
    run(rw, NAME)
    #rw.save(NAME)
    #data = load_results(NAME, **params)
    #print data.keys()
    
########### RETURN RW RESULTS, RUN IF NOT EXISTENT ###########
def get_results(NAME, params, calc=True):
    # check existing saved rws
    if fields.exists(NAME, **params):
        data = load_results(NAME, **params)
        N = len(data.times)
    else:
        N = 0
    # determine number of missing rws and run
    N_missing = params["N"] - N
    if N_missing > 0 and calc:
        new_params = Params(params, N=N_missing)
        rw = setup_rw(new_params)
        run(rw, NAME)
        rw.save(NAME)
    # return results
    data = load_results(NAME, **params)
    return data

########### PLOT EVOLUTION OF EXIT PROBABILITY ###########
def no_competitor_quantile(q):
    return -np.log(q)/275. * 1e9 # [ns]

q99 = no_competitor_quantile(.99)
q50 = no_competitor_quantile(.50)
endtime = 10e6
vlinecolor1 = "#999999"
vlinecolor2 = "#555555"
def plot_evolution(params, color=None, label=None,
                   vlines=False, vline_labels=False, xmax=endtime):
    data = get_results(NAME, params, calc=do_calculations)
    times = data.times
    success = data.success
    N = float(len(success))
    t = sorted(times[success])
    p = np.arange(sum(success))/N
    t.append(endtime)
    p = np.append(p, [p[-1]])
    errp = 1.96*np.sqrt(p*(1.-p)/N) # 95% confidence
    
    plt.semilogx(t, p, color=color, label=label)
    plt.fill_between(t, p - errp, p + errp, alpha=0.2,
                     facecolor=color, linewidth=0)
    
    plt.xlabel("Time [ns]")
    plt.ylabel("Exit probability")
    plt.xlim(xmin=0.1, xmax=xmax)
    
    if vlines:        
        plt.axvline(x=q99, linestyle=":", color=vlinecolor1, zorder=-100)
        plt.axvline(x=q50, linestyle="-", color=vlinecolor2, zorder=-100)
    
    print "last time: %.5f ms\nend prob: %.3f\nstd. dev.: %.3f" % (
        t[-2]*1e-6, p[-2], errp[-2])

def end_probability(params):
    data = get_results(NAME, params, calc=do_calculations)
    return data.success.mean()

if plot_rw_results:
    # FIGURE: Evolution for different starting positions
    Z = [0., 0.5, 1.0, 2.5, 5., 10.]
    P = []
    plt.figure("positions", figsize=(5, 4))
    for i, z in enumerate(Z):
        label = r"$z_0 = %.1f$nm" if z < 10 else r"$z_0 = %d$nm"
        plot_evolution(Params(params, zstart=z), xmax=1e6,
                       label=label % z, color="C%d" % i)
    plt.ylim(ymin=0.)
    plt.legend(frameon=False, loc="upper left")
    
    # FIGURE: Starting position vs. end probability
    zstop = params.zstop
    Z = [zstop, -1., -.75, -.5, -.25, 0., 0.25, 0.5, 1.0, 1.5, 2.5, 3.5, 5., 7.5, 10.]
    plt.figure("end_prob", figsize=(3, 4))
    P = [end_probability(Params(params, zstart=z)) for z in Z]
    plt.plot(Z, P, "o")
    # annotation of recognition site
    xofftext = -2.7
    yofftext = -0.15
    htext = 0.05
    color = "C1"
    plt.axvline(x=zstop, linestyle="--", color=color)
    plt.annotate("Recogni-\ntion site", (zstop, htext),
        xytext=(zstop + xofftext, htext + yofftext), color=color,
        )#arrowprops=dict(arrowstyle="->", color=color))
    #htext = 0.95
    #color = "#66aa66"
    #plt.axvline(x=0., linestyle="--", color=color)
    #plt.annotate("Pore entrance", (0., htext),
    #    xytext=(0. + xofftext, htext + yofftext), color=color,)
    #    #arrowprops=dict(arrowstyle="->", color=color))
    plt.xlabel(r"Distance $z_0$ [nm]")
    #plt.ylabel("Final exit probability")
    
    # FIGURE: Different binding durations
    plt.figure("bind_times", figsize=(4.5, 4))
    T = [0., 100., 10000.] #, 1000000.]
    labels = ["No binding", "100ns", u"10µs"] #, "1ms"]
    for i, t in enumerate(T):
        #label = r"$\tau = %g$ns" % t if t > 0 else "No binding"
        plot_evolution(Params(params, t_bind=t),
                       label=labels[i], color="C%d" % i,
                       vlines=(i == 0))
    htext = 0.25
    plt.text(q99 * 1.2, htext, "1%\n",  color=vlinecolor1)
    plt.text(q50 / 1.2, htext, "50%\n",  color=vlinecolor2, horizontalalignment='right')
    plt.text(np.sqrt(q50 * q99), htext, "\nconflicts",  color=vlinecolor1, horizontalalignment='center')
    #plt.xlim(xmin=1.)
    #plt.ylim(ymin=0.4)
    plt.ylim(ymin=0.)
    plt.legend(frameon=False, loc="upper left")#loc="lower right")
    
    # FIGURE: Different surface charges
    # TODO: surf charge together with bonding
    plt.figure("surf_charge", figsize=(4.5, 4))
    Rho = [-0.2, -0.1, 0.0001, 0.1, 0.2][::-1]
    labels = ["None"] # ["No binding", "100ns", u"10µs", "1ms"]
    for i, rho in enumerate(Rho):
        label = r"$\rho$ = %.1f C/m$^2$" % rho if rho is not None else "None"
        plot_evolution(Params(params, ahemqs=rho, ahemuniformqs=True), xmax=1e6,
                       label=label, color="C%d" % i)
    plt.ylim(ymin=0.)
    plt.ylabel("")
    plt.legend(frameon=False, loc="upper left")
    
    # FIGURE: surf charge vs. end prob, different bV
    plt.figure("surf_end_prob", figsize=(3, 4))
    V = [0.5, 1.0] #[1e-5, 0.25, 0.5, 1.0] # rho -0.3 + bV 1e-5, 0.1 doesnt converge
    Rho = [-0.2, -0.15, -0.1, -0.05, 0.0001, 0.025, 0.05, 0.1, 0.15, 0.2] # -0.3,
    colors = dict(zip(V, ["C2", "C3"]))
    for v in V:
        P = [end_probability(Params(params, bV=v, ahemqs=rho,
                                    ahemuniformqs=True)) for rho in Rho]
        print "Exit probs at %dmV" % (v*1000)
        print P
        print
        plt.plot(Rho, P, ":o", color=colors[v], label="%dmV" % (v*1000))
    plt.xlabel(r"Surface charge [C/m$^2$]")
    plt.legend(frameon=False)
    
    # FIGURE: Different applied voltages
    plt.figure("bV", figsize=(5, 4))
    V = [1e-5, 0.25, 0.5, 1.0]
    labels = ["None"] # ["No binding", "100ns", u"10µs", "1ms"]
    for i, v in enumerate(V):
        label = r"bV = %dmV" % (v*1000)
        plot_evolution(Params(params, bV=v), label=label, xmax=1e6,
                       color="C%d" % i)
    plt.ylim(ymin=0.)
    plt.legend(frameon=False, loc="upper left")
    
    # FIGURE: voltage vs. end prob
    plt.figure("bV_end_prob", figsize=(3, 4))
    #V = [-0.25, -0.1, 1e-5, 0.1, 0.2, 0.3, 0.4, 0.5, 1.]
    V = [-0.35, -0.25, -0.1, 1e-5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
    P = [end_probability(Params(params, bV=v)) for v in V]
    plt.plot(V, P, ":o")
    plt.xlabel(r"Voltage bias [V]")
    
    # FIGURE: starting position, different voltages:
    plt.figure("pos_bV_end_prob", figsize=(4, 4))
    zstop = params.zstop
    Z = [zstop, -.75, -.5, -.25, 0., 0.25, 0.5, 1.0, 1.5, 2., 3., 4., 5., 6., 7.5, 10., 12.5]
    rstart = lambda z: (1. if z>=0. else 0.1)
    #Z = [zstop, -.5, 0., 0.5, 1.0, 1.5, 2.5, 3.5, 5., 7.5, 10.]
    V = [1e-5, 0.25, 0.5, 1.0]
    for v in V:
        P = [end_probability(Params(params, zstart=z, bV=v, rstart=rstart(z))) for z in Z]
        plt.plot(Z, P, ":o", label="%dmV" % (v*1000))
    plt.legend(frameon=False)
    # annotation of recognition site
    xofftext = -2.7
    yofftext = -0.15
    htext = 0.0
    color = "#555555"
    plt.axvline(x=zstop, linestyle="--", color=color)
    plt.annotate("Recogni-\ntion site", (zstop, htext),
        xytext=(zstop + xofftext, htext + yofftext), color=color,
        )
    plt.xlabel(r"Distance $z_0$ [nm]")
    
########### CREATE CURRENT TRACE ###########
# linear interpolation:
def evaluate_interpolation(x, X, Y):
    """assuming data pairs X, Y = f(X) where X are sorted (monotonically increasing),
    find y = f(x) by interpolating between the nearest x values.
    x can be single number or array, return value will be of same type."""
    # get index of first element in X larger than x
    i = np.searchsorted(X, x)
    # now interpolate linearly between (X[i-1], Y[i-1]) and (X[i], Y[i])
    # x = X[i-1] + t*(X[i] - X[i-1])
    t = (x - X[i-1])/(X[i] - X[i-1])
    y = Y[i-1] + t*(Y[i] - Y[i-1])
    return y

def interpolation(X, Y, extend_to_infty=True):
    if extend_to_infty:
        X = [-np.infty] + X + [np.infty]
        Y = [Y[0]] + Y + [Y[-1]]
    X = np.array(X)
    Y = np.array(Y)
    return lambda x: evaluate_interpolation(x, X, Y)

def plot_current(i, rw, J, *plot_args, **plot_params):
    z = rw.positions[i][:, 2]
    t = rw.timetraces[i]
    plt.plot(t, J(z), *plot_args, **plot_params)
    
if create_current_trace:
    # create current profile
    rw = setup_rw(params)
    margin = 5.
    Htop = rw.ztop + 3*margin
    Hbot = -rw.zbot + 3*margin
    R = (Htop + Hbot)*.5
    sim_params0 = dict(Nmax=4e4, h=0.5, rDPore=.3)
    sim_params1 = dict(Nmax=2e4, h=0.5, rDPore=1.)
    sim_params = sim_params1
    sim_params = dict(get_active_params(params), R=R, Htop=Htop, Hbot=Hbot,
                      geop = dict(R=R, Htop=Htop, Hbot=Hbot), **sim_params)
    
    Z = np.linspace(rw.zbot - margin, rw.ztop + margin, 50)
    data = Iz(Z, nproc=5, name="current_exittime",
              calc=do_calculations, **sim_params)
    
    # FIGURE: current profile
    plt.figure("current")
    plt.plot(data.x, 0.3*np.array(data.J)*1e12, "-o")
    plt.xlabel("z position of nucleotide [nm]")
    plt.ylabel("Current [pA]")
    
    # define linear interpolation
    J = lambda x: 0.3e12*interpolation(data.x, data.J)(x)
    #zlin = np.linspace(rw.zbot - margin, rw.ztop + margin, 500)
    #plt.plot(zlin, 0.3e12*J(zlin), ".")
    
    # run a few random walks with recorded positions
    N = 20
    new_params = dict(params, N=N, record_positions=True, zstart=4.,
                      margtop=10, margbot=2)
    rw = get_rw("rw_exittime_path", new_params)
    zstop = params.zstop
    rzstop = 0.8
    # FIGURE(S): Molecule path + current trace
    stopcolor = "#ff6666"
    startcolor = "#77ee77"
    dotsize = 100
    max_num_figures = 5
    jsuccess = jfail = 0
    for i in range(N):
        x = rw.positions[i][:, 0]
        z = rw.positions[i][:, 2]
        t = rw.timetraces[i]
        if rw.success[i]:
            jsuccess += 1
            if jsuccess > max_num_figures: continue
            i0 = np.where(z < zstop)[0][0]
            plt.figure("current_trace_success%d" % jsuccess, figsize=(4, 4), frameon=False)
            rw.plot_path(i)
            plt.scatter([rw.positions[i][i0, 0]], [z[i0]], s=dotsize, c=stopcolor)
        else:
            jfail += 1
            if jfail > max_num_figures: continue
            plt.figure("current_trace_fail%d" % jfail, figsize=(4, 4), frameon=False)
            rw.plot_path(i)
        plt.xlim(-6, 17)
        plt.ylim(ymin=-16, ymax=10)
        plt.scatter([x[0]], [z[0]], c=startcolor, s=dotsize, label="Start")
        plt.plot([-rzstop, rzstop], [zstop, zstop], c=stopcolor,
                 label="Recognition site", linestyle="-", zorder=-50)
        print "Start", x[0], z[0]
        ax = plt.gca()
        ax.set_axis_off() # <=> plt.axis("off")
        ax.get_xaxis().set_visible(False) # so that white space disappears!
        ax.get_yaxis().set_visible(False)
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend(handles[::-1], labels[::-1], loc="lower left", frameon=False)
        #plt.gcf().tight_layout()
        #plt.subplots_adjust(left=0.5, right=0.95)
        ax = plt.axes([.57, .49, .27, .27])
        plot_current(i, rw, J, "-k")
        plt.ylim(300, 985)
        if rw.success[i]:
            plt.scatter([t[i0]], [J(z[i0])], s=dotsize, c=stopcolor)
        scalebar = AnchoredSizeBar(ax.transData, 10, "10ns", 2, 
                           pad=0.4,
                           color="#666666",
                           frameon=False,
                           size_vertical=25,
                           fontproperties = fm.FontProperties(size=8)
                           #fontproperties=fontprops
                           )
        ax.add_artist(scalebar)
        plt.xlabel("Time")
        plt.ylabel("Current")
        plt.tick_params(axis="both", which='both', left="off", labelleft="off",
                        bottom='off', labelbottom='off')
        #plt.gcf().tight_layout()
    
nanopores.savefigs("exittime", FIGDIR + "/ahem", ending=".pdf", bbox_inches="tight")
plt.show()
