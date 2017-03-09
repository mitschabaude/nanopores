# (c) 2017 Gregor Mitscha-Baude
from matplotlib import pyplot as plt
import numpy as np
import dolfin
from nanopores.tools import fields
fields.set_dir_dropbox()
from nanopores.models.nanopore import Setup
from nanopores.geometries.alphahempoly import poly
from nanopores.geometries.alphahem import default
from nanopores.geometries.cylpore import Pore, get_geo
from nanopores.models.diffusion_ahem import diff_profile_z_ahem, get_diffusivity

# params for precomputed diffusivity
params = dict(dim=2, Nmax=1e5, h=.5, ahemqsuniform=True, rMolecule=0.11)

#ap1 = 18
#ap2 = 49
#x0 = poly[18]
#x1 = poly[49]
#
#zmem = .5*(x0[1] + x1[1])
#print zmem
#
#poly = [[x[0], x[1] - zmem] for x in poly]
#proteincs = [z - zmem for z in default["proteincs"]]
#cs = [z - zmem for z in default["cs"]]
#default.update(zmem=0., hmem=2.82, Htop=10, Hbot=6, R=6, proteincs=proteincs, cs=cs)
#print default
#
#def new_get_geo(**params):
#    return get_geo(poly, **params)
#
#p = Pore(poly, **default)
#p.build(h=.5)
#
#p.polygons["alphahem"].plot("ok")
#p.polygons["membrane"].plot()
#p.polygons["bulkfluid_top"].plot()
#p.polygons["bulkfluid_bottom"].plot()
#plt.show()

#setup = Setup(get_geo=new_get_geo, geop=default, h=.5)
#setup = Setup(h=.5)
#setup.geo.plot_boundaries()
functions, mesh = fields.get_functions(name="Dalphahem-coupled", **params)
dist = functions["dist"]

#dolfin.plot(dist, interactive=True)

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

def diff_profile_fit(a=-10.3, b=0.05, N=20):
    Z = np.linspace(a, b, N)
    return Z, [Dfit(z) for z in Z]

z, D = diff_profile_fit(a=-12, b=2, N=100)
plt.plot(z, D, label="tabulated (infinite cylinder)")

data = diff_profile_z_ahem(**params)
z = [x0[2] for x0 in data["x"]]
Dz = data["D"]

plt.plot(z, Dz, "o-", label="computed with hydrodynamic model")
plt.ylabel("rel. diffusivity")
plt.xlabel("z [nm]")
plt.legend(loc="best")
#print results