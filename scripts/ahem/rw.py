# -*- coding: utf-8 -*-
# (c) 2017 Gregor Mitscha-Baude
"script to generate all figures for exit-time paper"
import numpy as np
import matplotlib.pyplot as plt
import nanopores
from nanopores.models.randomwalk import (get_pore, RandomWalk, run, load_results)
from nanopores.tools import fields
from nanopores import Params
fields.set_dir_dropbox()
FIGDIR = nanopores.dirnames.DROPBOX_FIGDIR

# TODO: Idee waere statt der ahem-validierung eine current trace abzubilden,
# damit man sieht an welchem punkt in der pore das Nukkleotid registriert wird
# oder eine echte current trace

########### PARAMETERS ###########
# TODO: confirm rMolecule and Qmol !!!
# goal: dt=0.01, T=10000, which means up to 1e6 steps
# take ~3.6s per 1000 samples, or ~10h for 100000 samples
params = nanopores.user_params(
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
do_calculations = True

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
    rw = setup_rw(params)
    run(rw, NAME, a=-11, b=-2, plot=True, save_count=1000)
    rw.save(NAME)
    
    data = load_results(NAME, **params)
    print data.keys()
    
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
        run(rw, NAME, plot=False)
        rw.save(NAME)
    # return results
    data = load_results(NAME, **params)
    return data

########### PLOT EVOLUTION OF EXIT PROBABILITY ###########
endtime = 5e6
def plot_evolution(params, color=None, label=None):
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
    plt.xlim(xmin=0.1, xmax=5e6)
    print "last time: %.5f ms\nend prob: %.3f\nstd. dev.: %.3f" % (
        t[-2]*1e-6, p[-2], errp[-2])

def end_probability(params):
    data = get_results(NAME, params, calc=do_calculations)
    return data.success.mean()

# FIGURE: Evolution for different starting positions
Z = [0., 0.5, 1.0, 2.5, 5., 10.]
P = []
plt.figure("positions", figsize=(5, 4))
for i, z in enumerate(Z):
    label = r"$z_0 = %.1f$nm" if z < 10 else r"$z_0 = %d$nm"
    plot_evolution(Params(params, zstart=z),
                   label=label % z, color="C%d" % i)
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
plt.figure("bind_times", figsize=(5, 4))
T = [0., 100., 10000., 1000000.]
labels = ["No binding", "100ns", u"10µs", "1ms"]
for i, t in enumerate(T):
    #label = r"$\tau = %g$ns" % t if t > 0 else "No binding"
    plot_evolution(Params(params, t_bind=t),
                   label=labels[i], color="C%d" % i)
#plt.xlim(xmin=1.)
#plt.ylim(ymin=0.4)
plt.legend(frameon=False, loc="upper left")#loc="lower right")

# FIGURE: Different surface charges
# TODO: surf charge together with bonding
plt.figure("surf_charge", figsize=(5, 4))
Rho = [-0.2, -0.1, 0.0001, 0.1, 0.2][::-1]
labels = ["None"] # ["No binding", "100ns", u"10µs", "1ms"]
for i, rho in enumerate(Rho):
    label = r"$\rho$ = %.1f C/m$^2$" % rho if rho is not None else "None"
    plot_evolution(Params(params, ahemqs=rho, ahemuniformqs=True),
                   label=label, color="C%d" % i)
plt.legend(frameon=False, loc="upper left")

# FIGURE: surf charge vs. end prob
plt.figure("surf_end_prob", figsize=(4, 4))
Rho = [-0.3, -0.2, -0.15, -0.1, -0.05, 0.0001, 0.025, 0.05, 0.1, 0.15, 0.2]
P = [end_probability(Params(params, ahemqs=rho, ahemuniformqs=True)) for rho in Rho]
plt.plot(Rho, P, "-o")
plt.xlabel(r"Surface charge [C/m$^2$]")

# FIGURE: Different applied voltages
plt.figure("bV", figsize=(5, 4))
V = [1e-5, 0.25, 0.5, 1.0]
labels = ["None"] # ["No binding", "100ns", u"10µs", "1ms"]
for i, v in enumerate(V):
    label = r"bV = %dmV" % (v*1000)
    plot_evolution(Params(params, bV=v), label=label, color="C%d" % i)
plt.legend(frameon=False, loc="upper left")

# FIGURE: voltage vs. end prob
plt.figure("bV_end_prob", figsize=(4, 4))
V = [-0.25, -0.1, 1e-5, 0.1, 0.2, 0.3, 0.4, 0.5, 1.]
P = [end_probability(Params(params, bV=v)) for v in V]
plt.plot(V, P, "-o")
plt.xlabel(r"Voltage bias [V]")

# FIGURE: starting position, different voltages:
plt.figure("pos_bV_end_prob", figsize=(3, 4))
zstop = params.zstop
Z = [zstop, -.5, 0., 0.5, 1.0, 1.5, 2.5, 3.5, 5., 7.5, 10.]
V = [1e-5, 0.25, 0.5, 1.0]
for v in V:
    P = [end_probability(Params(params, zstart=z, bV=v)) for z in Z]
    plt.plot(Z, P, "-o", label="bV = %dmV" % (v*1000))
plt.legend(frameon=False)
# annotation of recognition site
xofftext = -2.7
yofftext = -0.15
htext = 0.0
color = "black"
plt.axvline(x=zstop, linestyle="--", color=color)
plt.annotate("Recogni-\ntion site", (zstop, htext),
    xytext=(zstop + xofftext, htext + yofftext), color=color,
    )
plt.xlabel(r"Distance $z_0$ [nm]")

nanopores.savefigs("exittime", FIGDIR + "/ahem", ending=".pdf")
plt.show()
