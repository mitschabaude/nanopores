# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import folders
from nanopores.models.nanopore import IV
from matplotlib import pyplot as plt

ddata = dict(name="Dalphahem", dim=2, Nmax=.4e5, h=1., ahemqsuniform=True, rMolecule=0.11)
params = dict(dim=2, h=1., Nmax=2e4, rDPore=0.3, diffusivity_data=ddata)
V = [i/100. for i in range(-10, 11)]

#folders.fields.purge("IV-ahem")
results = IV(V, nproc=5, name="IV-ahem-pdD", ahemuniformqs=False, **params)
results_uniform = IV(V, nproc=5, name="IV-ahem-pdD", ahemuniformqs=True, **params)

plt.axvline(x=0, color="#999999", linestyle="--")
plt.axhline(y=0, color="#999999", linestyle="--")

V = 1e3*np.array(results["x"])
I = 1e12*np.array(results["J"])
plt.plot(V, I, "-", label="Simulation")
plt.xlabel("Voltage Bias [mV]")
plt.ylabel("Current [pA]")

V = 1e3*np.array(results_uniform["x"])
I = 1e12*np.array(results_uniform["J"])
plt.plot(V, I, "-", label="Simulation (homog. charge)")

# experimental data from Bhattacharya2011
bmV =  [-100., -71.42857143, -42.85714286, -14.28571429, 14.28571429, 42.85714286, 71.42857143, 100.]
I = [-83.339267043817102, -61.818625190008177, -39.496708569111611, -14.066625775593586, 14.6512949728476, 44.99789318249762, 76.122715987300211, 107.67609119161745]
plt.plot(bmV, I, "s", label="Experiment")
plt.legend(loc="best", frameon=False)

#plt.grid()

#plt.show()

from nanopores import savefigs
savefigs("ahemIV-pdD", folders.FIGDIR, size=(5,3.7))