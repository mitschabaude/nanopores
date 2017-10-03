# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker
import nanopores
import nanopores.models.pughpore as pugh
from nanopores.models.diffusion_interpolation import (cache_pugh_diffusivity,
                                                      diff_profile_z_pugh)

dparams = {2: dict(Nmax=1e5, dim=2, rMolecule=0.11, h=1.0),
           3: dict(Nmax=2e6, dim=3, rMolecule=0.11, h=2.0)}

# I for different pore diameters
@pugh.solvers.cache_forcefield("pugh_Idiam4_", pugh.defaultp)
def Idiam(diam, **params):
    dim = params["dim"]
    x0 = params.pop("x0")
    params.pop("diamPore")
    Jon = []
    Joff = []

    for d in diam:
        # get diffusivity interpolation
        cache_pugh_diffusivity(geoname="pugh2", diamPore=d, **dparams[dim])

        # current WITHOUT molecule
        setup = pugh.Setup(x0=None, diamPore=d, **params)
        pb, pnps = pugh.solve(setup, visualize=True)
        Jon.append(pnps.evaluate(setup.phys.CurrentPNPS)["J"])

        # current WITH molecule
        setup = pugh.Setup(x0=x0, diamPore=d, **params)
        pb, pnps = pugh.solve(setup, visualize=True)
        Joff.append(pnps.evaluate(setup.phys.CurrentPNPS)["J"])

    return dict(Jon=Jon, Joff=Joff)

params = {2: dict(dim=2, h=1., Nmax=1e5, x0=[0.,0.,0.], bV=-0.08,
                  diffusivity="Dpugh2"),
          3: dict(dim=3, h=2., Nmax=6e5, x0=[0.,0.,0.], bV=-0.08,
                  diffusivity="Dpugh2", stokesiter=True, cheapest=False)}

diam = {2: [4.18, 4.25, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.93, 5.2, 5.5, 6., 6.5, 7., 7.5],
        3: [4.22, 4.3, 4.4, 4.5, 4.6, 4.8, 5., 5.5, 6., 6.5, 7., 7.5]}
        #3: [4.18, 4.2, 4.25, 4.28, 4.3, 4.35, 4.4, 4.5, 4.555, 4.6, 4.8, 5., 5.5, 6., 6.5, 7., 7.5, 8.]}

fit = {2: 4.93, 3: 4.4}

calc = nanopores.user_param(calc=False)

if calc:
    for dim in 2, 3:
        for d in diam[dim]:
            diff_profile_z_pugh(nproc=6, diamPore=d)
#diam = {2: [3.7, 3.8, 3.9, 4.0, 4.141, 4.4, 4.6,
#            4.8, 5., 5.5, 6., 6.65],# 7.], # 7.5, 8.],
#        3: [4.18, 4.23, 4.3, 4.4,
#            4.5, 4.6, 4.8, 5.2, 5.5, 6., 7., 7.501]}

# semi-bad values: 4.17, 4.19, 4.21, 4.22, 4.25, 4.275, 4.7,
# bad values: 4.19, 4.9, 4.995, 5.095,

# for old setup:
#        3: [4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 4.23, 4.25, 4.275, 4.3, 4.4,
#            4.45, 4.5, 4.6, 4.8, 5., 5.5, 6., 6.5, 7., 7.5, 8.]}

figsize1 = (4.1, 3.3)
figsize2 = (4.1, 3.3)
for dim in 2, 3:
    fig = plt.figure("abs_%dD" % dim, figsize=figsize1)
    result = Idiam(diam[dim], calc=calc, nproc=2, **params[dim])
    d = result["x"]
    print "diameters (%dD):" % dim, d
    print "missing:", set(diam[dim]) - set(d)

    n = len(d)
    Jon = 1e12*np.array(result["Jon"])
    Joff = 1e12*np.array(result["Joff"])
    plt.plot(d, Jon, "s-b", label="No protein")
    plt.plot([0,10], [2.29*1e3*0.08]*2, "--b", label="Pugh et al.")
    plt.fill_between([0,10], [(2.29 - 0.26)*0.08*1e3]*2,
                     [(2.29 + 0.26)*0.08*1e3]*2, color="#ccccff")
    plt.plot(d, Joff, "s-g", label="With protein")
    plt.plot([0,10], [2.29*0.08*1e3*(1 - 0.262)]*2, "--g", label="Pugh et al.")
    plt.fill_between([0,10], [(2.29*(1 - 0.262) - 0.26)*0.08*1e3]*2,
                        [(2.29*(1 - 0.262) + 0.26)*0.08*1e3]*2, color="#ccffcc")
    plt.xlabel("Pore diameter [nm]")
    plt.ylabel("Current [pA]")

    plt.xlim(4., 7.6)
    plt.ylim(0., 1250.)
    locx, _ = plt.xticks()
    plt.axvline(x=2*2.0779, linestyle="--", color="#666666")
    plt.annotate("Diameter of trypsin", (2*2.0779, 900),
                 xytext=(2*2.0779 + 0.2, 900-20), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))
    #plt.title("influence of pore diameter on current (%dD)" %dim)
    plt.legend(loc="lower right")
    fig.tight_layout()

    fig = plt.figure("drop_%dD" % dim, figsize=figsize2)
    drop = (1. - Joff/Jon)*100
    plt.plot(d, drop, "s-r", label="Simulated blockade")
    plt.plot([0,10], [26.2]*2, "--r", label="Pugh et al.")
    plt.fill_between([0,10], [26.2 - 0.7]*2, [26.2 + 0.7]*2,
                     color="#ffcccc")
    plt.xlim(4., 7.6)
    #plt.xticks([4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5])
    #plt.ylim(ymin=0., ymax=50.)
    #plt.xlabel("pore diameter [nm]")
    plt.xlabel("Pore diameter [nm]")
    plt.ylabel("Current blockade [%]")
#    else:
#        loc, _ = plt.yticks()
#        plt.yticks(loc, [])
    ymax = 100 if dim==2 else 50
    htext = ymax*0.9
    plt.ylim(ymin=0., ymax=ymax)

    #loc, _ = plt.xticks()
    #plt.xticks(locx, [])
    #plt.yscale("log")
    #ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*0.01))
    #plt.gca().yaxis.set_major_formatter(ticks)

    plt.axvline(x=fit[dim], ymin=0., ymax=26.2/ymax, color="#ffaaaa", zorder=-90)
    plt.scatter([fit[dim]], [26.2], s=400, c="#ffaaaa", linewidths=0)
    plt.annotate("%.1f nm" % (fit[dim],), (fit[dim], 26.2),
             xytext=(fit[dim] + 0.1, 26.2 + 2.), color="#ff6666")

    plt.axvline(x=2*2.0779, linestyle="--", color="#666666")
    #htext = 45 # 2
    xofftext = 0.3 # 0.2
    plt.annotate("Diameter of trypsin", (2*2.0779, htext),
             xytext=(2*2.0779 + xofftext, htext-1.), color="#666666",
             arrowprops=dict(arrowstyle="->", color="#666666"))

#    plt.annotate("diameter of trypsin", (2*2.0779, 2),
#                 xytext=(2*2.0779 + 0.2, 2-0.2), color="#666666",
#                 arrowprops=dict(arrowstyle="->", color="#666666"))
#
    plt.legend(bbox_to_anchor=(.36, (.5 if dim==3 else .7)))
    fig.tight_layout()
    plt.subplots_adjust(left=0.2)
    #plt.legend(loc="center right" if dim==3 else "best")
    #fig.savefig(DIR + name + "_" + label + ".eps", bbox_inches=bbox_inches)


import folders
pugh.nano.savefigs("new", folders.FIGDIR + "/Idiam2", bbox_inches=None)
#plt.show()