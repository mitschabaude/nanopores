# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
import nanopores.models.pughpore as pugh
from nanopores.tools.utilities import collect_dict, user_param
import folders

# IV curve
@pugh.solvers.cache_forcefield("pugh_IV2", pugh.defaultp)
def IVDetail(V, **params):
    params["x0"] = None
    for bV, result in collect_dict(V):
        params["bV"] = bV
        setup = pugh.Setup(**params)
        pb, pnps = pugh.solve(setup)
        result.new = pnps.evaluate(setup.phys.CurrentPNPSDetail)
    return result

# I for different surface charges
@pugh.solvers.cache_forcefield("pugh_Irho2", pugh.defaultp)
def IrhoDetail(Rho, **params):
    params["x0"] = None
    for rho, result in collect_dict(Rho):
        params["dnaqsdamp"] = rho
        setup = pugh.Setup(**params)
        pb, pnps = pugh.solve(setup, visualize=True)
        result.new = pnps.evaluate(setup.phys.CurrentPNPSDetail)
    return result


ddata = {2: dict(name="Dpugh", Nmax=1e5, dim=2, r=0.11, h=1.0),
         3: dict(name="Dpugh", Nmax=2e6, dim=3, r=0.11, h=2.0)}

params = {2: dict(dim=2, h=1., Nmax=2e4, diffusivity_data=ddata[2]),
          3: dict(dim=3, h=2., Nmax=2e5, diffusivity_data=ddata[3],
                  stokesiter=False, cheapest=True)}

Rho = np.linspace(.1, .8, 8)
Rho = list(-Rho)[::-1] + [0.001] + list(Rho)

dnaqs = {3: -0.5882, 2: -0.6638}

def prep(J):
    return [-1e12*j for j in J[::-1]]

# TEST
#print IrhoDetail([-.5, .5], cache=False, dim=2, h=1., Nmax=1, diffusivity_data=ddata[2])

for dim in 2, 3:
    result = IrhoDetail(Rho, nproc=4, calc=user_param(calc=False), **params[dim])
    rho = result["x"]

    Jp = prep(result["Jp"])
    Jm = prep(result["Jm"])
    Jdif = prep(result["Jdif"])
    Jmig = prep(result["Jmig"])
    Jconv = prep(result["Jconv"])
    Jdm = np.array(Jdif) + np.array(Jmig)
    J = prep(result["J"])

    plt.figure("charge_%dD" % dim)
    plt.plot(rho, Jp, "s-", label=r"$J^+$", color="blue")
    plt.plot(rho, Jm, "s-", label=r"$J^-$", color="red")
    plt.xlabel("DNA surface charge [q/nm^2]")
    plt.ylabel("current  at -100 mV [pA]")
    plt.axvline(x=dnaqs[dim], linestyle="--", color="#666666")
    #plt.gca().yaxis.tick_right()
    plt.xlim(-.82, .82)
    plt.ylim(-200, 1100)
    plt.legend(loc="best")
    plt.annotate("estimated surface charge", (dnaqs[dim], 0),
                 xytext=(dnaqs[dim] + 0.1, 0-20), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))

    plt.figure("type_%dD" % dim)
    plt.plot(rho, Jmig, "s-", label=r"migrative", color="red")
    plt.plot(rho, Jconv, "o-", label=r"convective", color="blue")
    plt.plot(rho, Jdif, "v-", label=r"diffusive", color="green", zorder=-10)
    #plt.plot(rho, Jdm, "s-", label=r"diff. + mig.", color="cyan")
    #plt.plot(rho, J, "s-", label=r"total")
    plt.xlabel("DNA surface charge [q/nm^2]")
    plt.yticks(plt.yticks()[0], [])
    #plt.ylabel("current  at -100 mV [pA]")
    plt.axvline(x=dnaqs[dim], linestyle="--", color="#666666")
    plt.xlim(-.82, .82)
    plt.ylim(-200, 1100)
    plt.legend(loc="best")

pugh.nano.savefigs("Irho/detail", folders.FIGDIR, (5, 3.75))
#plt.show()