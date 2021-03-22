# (c) 2017 Gregor Mitscha-Baude
from matplotlib import pyplot as plt
import numpy as np
import dolfin
from nanopores.tools import fields
fields.set_dir_mega()
from nanopores.geometries.allpores import get_pore
from nanopores.models.nanopore import Setup
from nanopores.geometries.alphahempoly import poly
from nanopores.geometries.alphahem import default
from nanopores.geometries.cylpore import Pore, get_geo
from nanopores.models.diffusion_interpolation import Dt_plane
from nanopores.models.diffusion_ahem import diff_profile_z_ahem, get_diffusivity
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams.update({
    "font.size" : 7,
    "axes.titlesize" : 7,
    "font.family" : "sans-serif",
    "font.sans-serif" : ["CMU Sans Serif"],
    "lines.linewidth" : 1,
    "lines.markersize" : 2.5,
})

# params for precomputed diffusivity
params = dict(dim=2, Nmax=1e5, h=.5, ahemqsuniform=True, rMolecule=0.11)
functions, mesh = fields.get_functions(name="Dalphahem-coupled", **params)
dist = functions["dist"]

# construct D fit from Noskov2004 and plot tabulated D values
A = 0.64309
B = 0.00044
C = 0.06894
D = 0.35647
E = 0.19409

def Dfit(z, rion=0.11):
    rpore = dist([0., z])
    beta = rion/rpore
    return 1./(A + B*np.exp(beta/C) + D*np.exp(beta/E))

def Dr(z, rion=0.11):
    rpore = dist([0., z])
    return Dt_plane(rpore, rion)

z = np.linspace(-12, 2, 100)

poreRadii = [dist([0., zz]) for zz in z]
print("max", max(poreRadii))
print("min", min(poreRadii))

zbottom = get_pore(geoname="alphahem").protein.zmin()[1]
zNoskov = np.linspace(zbottom, 0, 100)

r = params['rMolecule']
Dr_ = [Dr(zz, r) for zz in z]
plt.plot(z, Dr_, "-b", label="r-dependent")

Dfit_ = [Dfit(zz) for zz in zNoskov]
plt.plot(zNoskov, Dfit_, "-", lw=1.5, color="orange", alpha=0.8, label="Noskov et al.")

data = diff_profile_z_ahem(a=-12, b=2, N=100, **params)
z1 = [x0[2] for x0 in data["x"]]
Dz = data["D"]
#plt.plot(z1, Dz, "-g", alpha=0.6)
plt.plot(z1, Dz, "og", label="LRNH", zorder=-5)

plt.ylabel("Rel. diffusivity")
plt.xlabel("Height of ion relative to pore top [nm]")
plt.xlim(-12, 2)
plt.ylim(0.5, 1.)
ax = plt.gca()
#ax.yaxis.tick_right()
#ax.yaxis.set_label_position("right")
plt.legend(loc="lower right", frameon=False)
plt.gcf().set_size_inches(3, 2.1)

from nanopores import savefigs
from folders import FIGDIR_HOWORKA
savefigs("Dz", FIGDIR_HOWORKA + "/ahem", ending=".pdf")
#print results
