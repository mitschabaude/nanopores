# (c) 2016 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
import nanopores.models.pughpore as pugh
from nanopores.models.diffusion_interpolation import cache_pugh_diffusivity

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
          3: dict(dim=3, h=2., Nmax=5e5, rDPore=0.95, x0=[0.,0.,0.], diamDNA=2.5,
                  bV=-0.08, stokesiter=True, cheapest=False)}

# old setup
params1 = dict(params)
params1[3] = dict(dim=3, h=2., Nmax=2e5, rDPore=0.95, x0=[0.,0.,0.], diamDNA=2.5,
                  bV=-0.08, stokesiter=False, cheapest=True)

diam = {2: [3.7, 3.8, 3.9, 4.0, 4.2, 4.4, 4.6,
            4.8, 5., 5.5, 6., 6.5, 7., 7.5, 8.],
        3: [4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 4.23, 4.25, 4.275, 4.3, 4.4,
            4.5, 4.6, 4.8, 5.]}
# for old setup:
#        3: [4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 4.23, 4.25, 4.275, 4.3, 4.4,
#            4.45, 4.5, 4.6, 4.8, 5., 5.5, 6., 6.5, 7., 7.5, 8.]}

for dim in 2, 3:
    plt.figure("abs_%dD" % dim)
    result = Idiam(diam[dim], calc=True, nproc=4, **params[dim])
    d = result["x"]
    print "diameters (%dD):" % dim, d
    print "missing:", set(diam[dim]) - set(d)
    n = len(d)
    Jon = 1e9*np.array(result["Jon"])
    Joff = 1e9*np.array(result["Joff"])
    plt.plot(d, Jon, "s-b", label="without molecule")
    plt.plot(d, [2.29*0.08]*n, "--b", label="Pugh et al.")
    plt.fill_between(d, [(2.29 - 0.26)*0.08]*n, [(2.29 + 0.26)*0.08]*n,
                     color="#ccccff")
    plt.plot(d, Joff, "s-g", label="with molecule")
    plt.plot(d, [2.29*0.08*(1 - 0.262)]*n, "--g", label="Pugh et al.")
    plt.fill_between(d, [(2.29*(1 - 0.262) - 0.26)*0.08]*n,
                        [(2.29*(1 - 0.262) + 0.26)*0.08]*n, color="#ccffcc")
    plt.xlabel("pore diameter [nm]")
    plt.ylabel("current at -80mV [nA]")
    #plt.title("influence of pore diameter on current (%dD)" %dim)
    plt.legend(loc="best")

    plt.figure("drop_%dD" % dim)
    drop = (1. - Joff/Jon)*100
    plt.plot(d, drop, "s-", label="% drop with molecule")
    plt.plot(d, [26.2]*n, "--k", label="Pugh et al.")
    plt.fill_between(d, [26.2 - 0.7]*n, [26.2 + 0.7]*n, color="#cccccc")
    plt.xlabel("pore diameter [nm]")
    plt.ylabel("current drop [%]")
    plt.legend(loc="best")

import folders
pugh.nano.savefigs("pugh_Idiam", folders.FIGDIR, (10, 8))
#plt.show()