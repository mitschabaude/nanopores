# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import folders
fields = folders.fields
from nanopores.models.nanopore import IV
import nanopores.plots as plots
colors = plots.colors
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
adam1 = '#3069ab'
adam2 = '#5993d0'
adam3 = '#93cbf2'
intense = '#0224bd'

blue = colors.pure
#pink = "#9900ff" #plots.pure(colors.pink)

colors_list_hue = ['k',
               plots.set_sl("#0000ff", .7, .5),
               plots.set_sl("#00aaff", .8, .6),
               plots.set_sl("#00ffaa", .9, .8),]

colors_list = [plots.set_sl("#6600ff", .8, .7),
               plots.set_sl(blue, 1., .15),
               plots.set_sl(blue, 1., .6),
               plots.set_sl(blue, 1., .85),]

line_width = 1.5 # 1 normal, 1.5 thicker

ddsimplerough = dict(name="Dalphahem", dim=2, Nmax=.1e5, h=1.,
             ahemqsuniform=True, rMolecule=0.11)

ddsimplefine = dict(name="Dalphahem", dim=2, Nmax=1e5, h=.5,
             ahemqsuniform=True, rMolecule=0.11)

ddcoupled = dict(name="Dalphahem-coupled", dim=2, Nmax=1e5, h=.5,
             ahemqsuniform=True, rMolecule=0.11)

ddprofile = dict(name="Dalphahem-profile", dim=2, Nmax=1e5, h=.5,
             ahemqsuniform=True, rMolecule=0.11)

def get_IV(calc=True, **params):
    V = [i/100. for i in range(-10, 11)]

    #folders.fields.purge("IV-ahem")
    results = IV(V, nproc=5, name="IV-ahem", calc=calc, ahemuniformqs=False, **params)
    results_uniform = IV(V, nproc=5, name="IV-ahem", calc=calc, ahemuniformqs=True, **params)

    plot_grid()

    V = 1e3*np.array(results["x"])
    I = 1e12*np.array(results["J"])
    plt.plot(V, I, "-b", label="Simulation")
    plt.xlabel("Voltage Bias [mV]")
    plt.ylabel("Current [pA]")

    V = 1e3*np.array(results_uniform["x"])
    I = 1e12*np.array(results_uniform["J"])
    plt.plot(V, I, ":g", label="Simulation (homog. charge)")

    # experimental data from Bhattacharya2011
    bmV =  [-100., -71.42857143, -42.85714286, -14.28571429, 14.28571429, 42.85714286, 71.42857143, 100.]
    I = [-83.339267043817102, -61.818625190008177, -39.496708569111611, -14.066625775593586, 14.6512949728476, 44.99789318249762, 76.122715987300211, 107.67609119161745]
    plt.plot(bmV, I, "sr", label="Experiment", markersize=3)
    ax = plt.gca()
    ax.get_yaxis().set_ticks_position('right')
    ax.get_yaxis().set_label_position('right')
    plt.legend(loc="best", frameon=False)

def plot_grid(color="#aaaaaa"):
    plt.axvline(x=0, color=color, linestyle="-")
    plt.axhline(y=0, color=color, linestyle="-")

def plot_experiment():
    bmV =  [-100., -71.42857143, -42.85714286, -14.28571429, 14.28571429, 42.85714286, 71.42857143, 100.]
    I = [-83.339267043817102, -61.818625190008177, -39.496708569111611, -14.066625775593586, 14.6512949728476, 44.99789318249762, 76.122715987300211, 107.67609119161745]
    plt.plot(bmV, I, "s", color=colors.experiment, label="Experiment")
    #G = 1e-3*I[2]/(-0.04285) # -40
    G = 1e-3*I[5]/(0.04285) # +40
    print "Conductivity experimental: %.4f nS" % (G,)
    return G

def plot_experiment_simple():
    bmV =  [-100., -71.42857143, -42.85714286, -14.28571429, 14.28571429, 42.85714286, 71.42857143, 100.]
    I = [-83.339267043817102, -61.818625190008177, -39.496708569111611, -14.066625775593586, 14.6512949728476, 44.99789318249762, 76.122715987300211, 107.67609119161745]
    plt.plot(bmV, I, "s", color=colors.experiment, label="Experiment")
    #G = 1e-3*I[2]/(-0.04285) # -40
    G = 1e-3*I[5]/(0.04285) # +40
    return G

def fine_vs_coarse(calc=True, **params):
    V = [i/100. for i in range(-10, 11)]
    results_fine = IV(V, nproc=5, name="IV-ahem", calc=calc,
                      ahemuniformqs=False, diffusivity_data=ddsimplefine,
                      **params)
    results_rough = IV(V, nproc=5, name="IV-ahem", calc=calc,
                      ahemuniformqs=False, diffusivity_data=ddsimplerough,
                      **params)
    I_fine = 1e12*np.array(results_fine["J"])
    I_rough = 1e12*np.array(results_rough["J"])
    V = 1e3*np.array(results_fine["x"])
    plt.plot(V, I_fine, label="fine mesh for D")
    plt.plot(V, I_rough, label="coarse mesh for D")
    plt.xlabel("Voltage Bias [mV]")
    plt.ylabel("Current [pA]")
    plt.legend(loc="best", frameon=False)

def compare_D_models(calc=True, **params):
    params["ahemuniformqs"] = False
    V = [i/100. for i in range(-10, 11)]
    Vplot = 1e3*np.array(V)
    DD = OrderedDict([
        ("Bulk diff.", None),
        ("r-dependent", ddsimplefine),
        ("z-dependent", ddprofile),
        ("r- and z-dep.", ddcoupled),
    ])
    #colors = ["k", "b", "g", "c"]
    #colors_ = ["k", colors.pink, colors.darkintense, colors.muted]
    #dashes = [[1000,1], [6,2], [6,1,1,1], [6,1,1,1,1,1]]
    #lines = ["--", "-", ":", "-."]
    lines = ["-", "-", "-", "-"]
    plot_grid()
    G = [0]*len(DD)
    for i, model in enumerate(DD):
        params["diffusivity_data"] = DD[model]
        mod_params = dict(params)
        if model == "Bulk diff.":
            mod_params["rDPore"] = 1.
        results = IV(V, nproc=3, name="IV-ahem", calc=calc, **mod_params)
        I = 1e12*np.array(results["J"])
        plt.plot(Vplot, I, "-", color=colors_list[i], linestyle=lines[i],
                 label=model, lw=line_width)

        # print conductivities
        #G[i] = 1e-3*I[V.index(-0.04)]/(-0.04) # -40
        G[i] = 1e-3*I[V.index(0.04)]/(0.04) # +40
        print "Conductivity %s: %.4f nS" % (model, G[i])

    plt.xlabel("Voltage [mV]")
    plt.ylabel("Current [pA]")
    plt.ylim(-100, 100)
    plt.ylim(-200, 200)
    plt.xticks([-100, 0, 100])
    plt.yticks([-200, 0, 200])
    # plots.addMinorTicks()
    plots.removeTopRightFrame()
    gexp = plot_experiment()
    plt.legend(loc="upper left", frameon=False)

    x = 10
    y = -35
    plt.text(x, y, "Overest. (+40 mV):")
    #y -= 32
    #plt.text(x, y, "")
    for i, g in enumerate(G):
        y -= 36
        change = int(100*(g/gexp - 1.))
        plt.text(x, y, "+%d%%" % change, color=colors_list[i])
        
def compare_D_models_simple(calc=True, **params):
    params["ahemuniformqs"] = False
    V = [i/100. for i in range(-10, 11)]
    Vplot = 1e3*np.array(V)
    DD = OrderedDict([
        ("Simulation", ddcoupled),
    ])
    #colors = ["#00cc00"]
    #plot_grid()
    G = [0]*len(DD)
    for i, model in enumerate(DD):
        params["diffusivity_data"] = DD[model]
        mod_params = dict(params)
        if model == "Bulk D":
            mod_params["rDPore"] = 1.
        results = IV(V, nproc=3, name="IV-ahem", calc=calc, **mod_params)
        I = 1e12*np.array(results["J"])
        plt.plot(Vplot, I, "-", color=colors.mediumlight, label=model,
                lw=1.5, zorder=-200)

    plt.xlabel("Voltage [mV]")
    plt.ylabel("Current [pA]")
    plt.xlim(-115, 125)
    ax = plt.gca()
    ax.set_yticks([-100, 0, 100])
    ax.tick_params(axis="y", direction='in', length=2, pad=-5)
    ax.set_yticklabels(["", 0, 100], fontdict={"ha": "left"})
    # plt.tick_params(
    #     axis="both",          # changes apply to the x-axis
    #     which="both",      # both major and minor ticks are affected
    #     bottom=False,      # ticks along the bottom edge are off
    #     top=False,         # ticks along the top edge are off
    #     left=False,
    #     right=False,
    #     labelleft=False,
    #     labelbottom=False)
    #plt.ylim(-100, 100)
    #plt.ylim(-200, 200)
    plots.removeTopRightFrame(ax)
    plot_experiment_simple()
    
    plt.legend(loc="lower right", frameon=False, bbox_to_anchor=(1.07, 0))
    plt.title("IV curve")

calc = False
default = dict(geoname="alphahem", dim=2, h=1., Nmax=2e4, rDPore=0.3, ahemqs=-0.1)

# with constant D = 0.3*D0 in pore
dd = None
params = dict(default, diffusivity_data=dd)
plt.figure("rD03", figsize=(2.6, 2.))
get_IV(calc, **params)

# with 'simple' D interpolation (only analytical near-plane model)
dd = ddsimplefine
params = dict(default, diffusivity_data=dd)
plt.figure("pdDsimple")
get_IV(calc, **params)

# with analytical r-dependence of D multiplied by z profile
dd = ddcoupled
params = dict(default, diffusivity_data=dd)
plt.figure("pdDcoupled")
get_IV(calc, **params)

# comparison
#plt.figure("finevscoarse")
#fine_vs_coarse(False, **default)
plt.figure("compareD", figsize=(2.7, 2.1)) #, dpi=300)
compare_D_models(calc, **default)

plt.figure("compareD_simple", figsize=(1.7, 1.6))
compare_D_models_simple(calc, **default)

from nanopores import savefigs
#savefigs("IV", folders.FIGDIR_HOWORKA + "/ahem", ending=".pdf")
savefigs("IV", "./tmp/figures" + "/ahem", ending=".pdf")
#plt.show()
