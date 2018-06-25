# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import folders
fields = folders.fields
from nanopores.models.nanopore import IV
from collections import OrderedDict
from matplotlib import pyplot as plt

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
    plt.plot(V, I, "-g", label="Simulation (homog. charge)")

    # experimental data from Bhattacharya2011
    bmV =  [-100., -71.42857143, -42.85714286, -14.28571429, 14.28571429, 42.85714286, 71.42857143, 100.]
    I = [-83.339267043817102, -61.818625190008177, -39.496708569111611, -14.066625775593586, 14.6512949728476, 44.99789318249762, 76.122715987300211, 107.67609119161745]
    plt.plot(bmV, I, "sr", label="Experiment")
    plt.legend(loc="best", frameon=False)

def plot_grid():
    plt.axvline(x=0, color="#999999", linestyle="-")
    plt.axhline(y=0, color="#999999", linestyle="-")

def plot_experiment():
    bmV =  [-100., -71.42857143, -42.85714286, -14.28571429, 14.28571429, 42.85714286, 71.42857143, 100.]
    I = [-83.339267043817102, -61.818625190008177, -39.496708569111611, -14.066625775593586, 14.6512949728476, 44.99789318249762, 76.122715987300211, 107.67609119161745]
    plt.plot(bmV, I, "sr", label="Experiment")
    #G = 1e-3*I[2]/(-0.04285) # -40
    G = 1e-3*I[5]/(0.04285) # +40
    print "Conductivity experimental: %.4f nS" % (G,)
    return G

def plot_experiment_simple():
    bmV =  [-100., -71.42857143, -42.85714286, -14.28571429, 14.28571429, 42.85714286, 71.42857143, 100.]
    I = [-83.339267043817102, -61.818625190008177, -39.496708569111611, -14.066625775593586, 14.6512949728476, 44.99789318249762, 76.122715987300211, 107.67609119161745]
    plt.plot(bmV, I, "sr", markersize=8, label="Experiment")
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
        ("Bulk D", None),
        ("r-dependent D", ddsimplefine),
        ("z-dependent D", ddprofile),
        ("r- and z-dep. D", ddcoupled),
    ])
    colors = ["k", "b", "g", "c"]
    plot_grid()
    G = [0]*len(DD)
    for i, model in enumerate(DD):
        params["diffusivity_data"] = DD[model]
        mod_params = dict(params)
        if model == "Bulk D":
            mod_params["rDPore"] = 1.
        results = IV(V, nproc=3, name="IV-ahem", calc=calc, **mod_params)
        I = 1e12*np.array(results["J"])
        plt.plot(Vplot, I, "-", color=colors[i], label=model)

        # print conductivities
        #G[i] = 1e-3*I[V.index(-0.04)]/(-0.04) # -40
        G[i] = 1e-3*I[V.index(0.04)]/(0.04) # +40
        print "Conductivity %s: %.4f nS" % (model, G[i])

    plt.xlabel("Voltage [mV]")
    plt.ylabel("Current [pA]")
    plt.ylim(-100, 100)
    plt.ylim(-200, 200)
    gexp = plot_experiment()
    plt.legend(loc="best", frameon=False)

    y = -35
    plt.text(15, y, "Conductance ")
    y -= 30
    plt.text(15, y, "overestimate (+40 mV):")
    for i, g in enumerate(G):
        y -= 30
        change = int(100*(g/gexp - 1.))
        plt.text(15, y, "+%d%%" % change, color=colors[i])
        
def compare_D_models_simple(calc=True, **params):
    params["ahemuniformqs"] = False
    V = [i/100. for i in range(-10, 11)]
    Vplot = 1e3*np.array(V)
    DD = OrderedDict([
        ("Simulation", ddcoupled),
    ])
    colors = ["C0"]
    plot_grid()
    G = [0]*len(DD)
    for i, model in enumerate(DD):
        params["diffusivity_data"] = DD[model]
        mod_params = dict(params)
        if model == "Bulk D":
            mod_params["rDPore"] = 1.
        results = IV(V, nproc=3, name="IV-ahem", calc=calc, **mod_params)
        I = 1e12*np.array(results["J"])
        plt.plot(Vplot, I, "-", color=colors[i], label=model)

    plt.xlabel("Voltage [mV]")
    plt.ylabel("Current [pA]")
    #plt.ylim(-100, 100)
    #plt.ylim(-200, 200)
    gexp = plot_experiment_simple()
    plt.legend(loc="best", frameon=False)
    plt.title("IV curve")

calc = False
default = dict(geoname="alphahem", dim=2, h=1., Nmax=2e4, rDPore=0.3, ahemqs=-0.1)

# with constant D = 0.3*D0 in pore
dd = None
params = dict(default, diffusivity_data=dd)
plt.figure("rD03")
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
plt.figure("compareD", figsize=(5,3.7))
compare_D_models(calc, **default)

plt.figure("compareD_simple", figsize=(2.7, 2.3))
compare_D_models_simple(calc, **default)

from nanopores import savefigs
savefigs("IV", folders.FIGDIR + "/ahem")
#plt.show()