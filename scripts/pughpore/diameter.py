# (c) 2016 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import nanopores
import nanopores.models.pughpore as pugh
from nanopores.models.diffusion_interpolation import cache_pugh_diffusivity
from math import sqrt, pi

dparams = {2: dict(Nmax=1e5, dim=2, r=0.11, h=1.0),
           3: dict(Nmax=2e6, dim=3, r=0.11, h=2.0)}

# I for different pore diameters
@pugh.solvers.cache_forcefield("pugh_Idiam", pugh.defaultp)
def Idiam(diam, **params):
    dim = params["dim"]
    x0 = params.pop("x0")
    params.pop("diamPore")
    Jon = []
    Joff = []
    # TODO bad hack
    if dim==3:
        params.update(h=1.5, Nmax=7e5)
    if dim==2:
        params.update(h=0.75, Nmax=1e5)

    for d in diam:
        # get diffusivity interpolation
        cache_pugh_diffusivity(diamPore=d, **dparams[dim])
        ddata = dict(dparams[dim], name="Dpugh", diamPore=d)

        # current WITHOUT molecule
        setup = pugh.Setup(x0=None, diamPore=d, diffusivity_data=ddata, **params)
        pb, pnps = pugh.solve(setup, visualize=True)
        Jon.append(pnps.evaluate(setup.phys.CurrentPNPS)["J"])

        # current WITH molecule
        setup = pugh.Setup(x0=x0, diamPore=d, diffusivity_data=ddata, **params)
        pb, pnps = pugh.solve(setup, visualize=True)
        Joff.append(pnps.evaluate(setup.phys.CurrentPNPS)["J"])

    return dict(Jon=Jon, Joff=Joff)

params = {2: dict(dim=2, h=1., Nmax=2e4, rDPore=0.95, x0=[0.,0.,0.], diamDNA=2.5,
                  bV=-0.08),
          3: dict(dim=3, h=2., Nmax=6e5, rDPore=0.95, x0=[0.,0.,0.], diamDNA=2.5,
                  bV=-0.08, stokesiter=True, cheapest=False)}

# old setup
params1 = dict(params)
params1[3] = dict(dim=3, h=2., Nmax=2e5, rDPore=0.95, x0=[0.,0.,0.], diamDNA=2.5,
                  bV=-0.08, stokesiter=False, cheapest=True)

diam = {2: [3.7, 3.8, 3.9, 4.0, 4.141, 4.4, 4.6,
            4.8, 5., 5.5, 6., 6.65],# 7.], # 7.5, 8.],
        3: [4.18, 4.23, 4.3, 4.4,
            4.5, 4.6, 4.8, 5.2, 5.5, 6., 7., 7.501]}

# semi-bad values: 4.17, 4.19, 4.21, 4.22, 4.25, 4.275, 4.7,
# bad values: 4.19, 4.9, 4.995, 5.095,

# for old setup:
#        3: [4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 4.23, 4.25, 4.275, 4.3, 4.4,
#            4.45, 4.5, 4.6, 4.8, 5., 5.5, 6., 6.5, 7., 7.5, 8.]}
calc = nanopores.user_param(calc=False)

for dim in 2, 3:
    plt.figure("abs_%dD" % dim)
    result = Idiam(diam[dim], calc=calc,
                   nproc=4, **params[dim])
    d = result["x"]
    print "diameters (%dD):" % dim, d
    print "missing:", set(diam[dim]) - set(d)
    if dim==2:
        d = [dd*2./sqrt(pi) for dd in d]
    n = len(d)
    Jon = 1e12*np.array(result["Jon"])
    Joff = 1e12*np.array(result["Joff"])
    plt.plot(d, Jon, "s-b", label="without molecule")
    plt.plot([0,10], [2.29*1e3*0.08]*2, "--b", label="Pugh et al.")
    plt.fill_between([0,10], [(2.29 - 0.26)*0.08*1e3]*2,
                     [(2.29 + 0.26)*0.08*1e3]*2, color="#ccccff")
    plt.plot(d, Joff, "s-g", label="with molecule")
    plt.plot([0,10], [2.29*0.08*1e3*(1 - 0.262)]*2, "--g", label="Pugh et al.")
    plt.fill_between([0,10], [(2.29*(1 - 0.262) - 0.26)*0.08*1e3]*2,
                        [(2.29*(1 - 0.262) + 0.26)*0.08*1e3]*2, color="#ccffcc")
    plt.xlabel("pore diameter [nm]")
    plt.ylabel("current at -80mV [pA]")
    plt.xlim(4., 7.6)
    plt.ylim(0., 1250.)
    plt.axvline(x=2*2.0779, linestyle="--", color="#666666")
    plt.annotate("diameter of trypsin", (2*2.0779, 900),
                 xytext=(2*2.0779 + 0.2, 900-20), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))
    #plt.title("influence of pore diameter on current (%dD)" %dim)
    plt.legend(loc="lower right")

    plt.figure("drop_%dD" % dim)
    drop = (1. - Joff/Jon)*100
    plt.plot(d, drop, "s-r", label="% drop with molecule")
    plt.plot([0,10], [26.2]*2, "--r", label="Pugh et al.")
    plt.fill_between([0,10], [26.2 - 0.7]*2, [26.2 + 0.7]*2,
                     color="#ffcccc")
    plt.xlabel("pore diameter [nm]")
    plt.ylabel("current drop [%]")
    plt.xlim(4., 7.6)
    plt.ylim(ymin=0., ymax=100.)
    #plt.yscale("log")
    #ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*0.01))
    #plt.gca().yaxis.set_major_formatter(ticks)
    if dim==2:
        plt.axvline(x=4.673, ymin=0., ymax=26.2/100., color="#ffaaaa", zorder=-90)
        plt.scatter([4.673], [26.2], s=400, c="#ffaaaa", linewidths=0)
        plt.annotate("4.7 nm", (4.673, 26.2),
                 xytext=(4.673 + 0.1, 26.2 + 2.), color="#ff6666")
                 #arrowprops=dict(arrowstyle="->", color="#666666"))

    plt.axvline(x=2*2.0779, linestyle="--", color="#666666")
    htext = 80 # 2
    xofftext = 0.3 # 0.2
    plt.annotate("diameter of trypsin", (2*2.0779, htext),
                 xytext=(2*2.0779 + xofftext, htext-1.), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))

#    plt.annotate("diameter of trypsin", (2*2.0779, 2),
#                 xytext=(2*2.0779 + 0.2, 2-0.2), color="#666666",
#                 arrowprops=dict(arrowstyle="->", color="#666666"))
    plt.legend(loc="center right")
    #plt.legend(loc="center right" if dim==3 else "best")


import folders
pugh.nano.savefigs("pugh_Idiam", folders.FIGDIR, (6, 4.5))
#plt.show()