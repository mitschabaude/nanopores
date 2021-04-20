# (c) 2019 Gregor Mitscha-Baude
import folders
import numpy as np
import nanopores as nano
from nanopores.geometries import get_pore
from nanopores.models.diffusion import diffusivity_tensor
from nanopores.models.diffusion_interpolation import Dt_plane
from nanopores.tools import fields, collect_dict, plot_sliced
from nanopores.tools.solvers import cache_forcefield
from nanopores.models import nanopore
import nanopores.plots as plots
from matplotlib import pyplot as plt
fields.set_dir_mega()

showSolid = False

if showSolid:
    setup = nanopore.Setup(geoname="alphahem", dim=3, Nmax=5e4, h=1., geop=dict(subs='solid'))
    setup.geo.plot_subdomains()
    plot_sliced(setup.geo)

r = 0.11
params = dict(Nmax=2e5, h=.5, geop=dict(Htop=5., Hbot=13., R=8.), 
              z0=-7.5, rMolecule = r)

# evaluation points
def evaluation_points(N, **params):
    z = params['z0']
    pore = get_pore(geoname="alphahem", **params)
    a = pore.radius_at(z) - params['rMolecule']*1.5
    return [[x, 0., z] for x in np.linspace(a, 0, N)]

@cache_forcefield("diff_r_ahem")
def diff_r_ahem(X, **params):
    print 'POINTS: ' + repr(X)
    for x, result in collect_dict(X):
        print 'NEXT POINT: ' + repr(x)
        params["x0"] = x
        setup = nanopore.Setup(geoname="alphahem", dim=3, **params)
        D = diffusivity_tensor(setup)
        plot_sliced(setup.geo)
        result.new = dict(D=D)
    return result

# phys = nano.Physics("pore_mol")
# D00 = phys.kT/(6.*phys.pi*phys.eta*r*1e-9)
# print(D00)

X = evaluation_points(10, **params)
data = diff_r_ahem(X, **params);

R = get_pore(geoname="alphahem", **params).radius_at(params['z0'])
x = [R - xi[0] for xi in data["x"]]
DD = data['D']

Dxx = [D.load()[0][0] for D in DD]
Dyy = [D.load()[1][1] for D in DD]
Dzz = [D.load()[2][2] for D in DD]

#plt.plot(x, Dxx, "ob", label=r"$D_{xx}$")
#plt.plot(x, Dyy, "sg", label=r"$D_{yy}$")

Dxx1 = Dxx[-1]
Dyy1 = Dyy[-1]
Dzz1 = Dzz[-1]

xlin = np.linspace(r+1e-3, R, 100)
#dn = [Dn_plane(t, r, N=20) for t in xlin]
#dn = [d*Dxx1/dn[-1] for d in dn]
#plt.plot(xlin, dn, "-b")
#dt = [Dt_plane(t, r) for t in xlin]
#dt = [d*Dyy1/dt[-1] for d in dt]
#plt.plot(xlin, dt, "-g")

D0 = [1. for t in xlin]
D_r = [Dt_plane(t, r) for t in xlin]
D_z = [Dzz1 for t in xlin]
D_rz = [d*Dzz1/D_r[-1] for d in D_r]

MD = [(1. - np.exp(-(t - 0.22)/0.162)) for t in xlin]

plt.figure('main')
lw = 1.5
#plt.plot(xlin, MD, "-", label=r"r-dependent (MD)", color=plots.rotateRed(plots.colors.muted))
plt.plot(xlin, D_r, "-", label="r-dependent", color=plots.colors.dark, lw=lw)
plt.plot(xlin, D_z, "-", label="z-dependent", color=plots.colors.mediumintense, lw=lw)
plt.plot(xlin, D_rz, "-", label=r"r- and z-dep.", color=plots.colors.light, lw=lw)
plt.plot(x, Dzz, "o", label=r"LRNH", color=plots.colors.mediumintense)

plt.xlim(0., R+0.05)
plt.xticks([0, 0.5, 1.])

plt.xlabel("Ion distance from pore wall [nm]")
plt.ylabel("Rel. diffusivity")
plt.ylim(0, 1)

plt.axvline(x=0.11, linestyle="--", color="#666666")
plt.annotate("Ion radius", (0.11, 0.94),
                 xytext=(0.25, 0.94-0.002), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))
#plt.yticks([i/10. for i in range(0, 11, 2)])
plt.yticks([0, .5, 1])

plots.removeTopRightFrame()
plt.legend(loc="lower right", frameon=False)
plt.gcf().set_size_inches(2.7, 2.1)

plt.figure('MD')
plt.plot(xlin, MD, "-", label=r"r-dependent (MD)", color=plots.rotateRed(plots.colors.muted))
plt.plot(xlin, D_r, "-", label="r-dependent (LRNH)", color=plots.colors.dark)

#plt.plot(xlin, D_z, "-", label="z-dependent", color=plots.colors.medium)
#plt.plot(xlin, D_rz, "-", label=r"r- and z-dep.", color=plots.colors.lightmuted)
#plt.plot(x, Dzz, "o", label=r"LRNH", color=plots.colors.medium)

plt.xlim(0., R+0.05)
plt.xticks([0, 0.5, 1.])

plt.xlabel("Ion distance from plane wall [nm]")
plt.ylabel("Rel. diffusivity")
plt.ylim(0, 1)

plt.axvline(x=0.11, linestyle="--", color="#666666")
plt.annotate("Ion radius", (0.11, 0.94),
                 xytext=(0.25, 0.94-0.002), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))
#plt.yticks([i/10. for i in range(0, 11, 2)])
plt.yticks([0, .5, 1])

plots.removeTopRightFrame()
plt.legend(loc="lower right", frameon=False)
plt.gcf().set_size_inches(2.7, 2.1)

from nanopores import savefigs
savefigs("Dr", folders.FIGDIR_CWD + "/ahem", ending=".pdf")
