# (c) 2018 Gregor Mitscha-Baude
"plot simulated pores schematically with protein/receptors"
from matplotlib import rcParams, rc
rcParams.update({
    "font.size" : 7,
    "axes.titlesize" : 7,
    "font.family" : "sans-serif",
    "font.sans-serif" : ["CMU Sans Serif"],
    "lines.linewidth" : 1,
    "lines.markersize" : 5,
})
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from nanopores import get_pore, savefigs
from nanopores.models.randomwalk import RandomWalk, Ball

cprot = "#bd788a"
crec = "#788abd"

# for plotting proteins andreceptors
ball_settings = dict(facecolor="#eeeeee", edgecolor="k", linewidth=1.)
def add_ball(ax, x, y, r, **settings):
    p = mpatches.Circle((x, y), r, **dict(ball_settings, **settings))
    p.set_zorder(200)
    ax.add_patch(p)

###############################################################

plt.figure("ahem", figsize=(1.7, 1.7))
pore = get_pore(geoname="alphahem", R=8, Htop=3.5, Hbot=12.5)
rw = RandomWalk(pore, rMolecule=0.11)
ax = rw.plot_pore()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.axis('off')
plt.title("Protein pore")

###############################################################

plt.figure("wei", figsize=(1.7, 1.7))
pore = get_pore(geoname="wei",
    R=100, Htop=120, Hbot=80)
rw = RandomWalk(pore, rMolecule=3)

rrec = 5
zrec = rw.zbot + rrec + (rw.ztop - rw.zbot - 2.*rrec)*0.5
xrec = pore.radius_at(zrec) - rrec - 1.

rprot = 10
zprot = rw.ztop + 20.
xprot = -10.

ax = rw.plot_pore()
add_ball(ax, xrec, zrec, rrec, facecolor=crec, linewidth=0)
add_ball(ax, xprot, zprot, rprot, facecolor=cprot, linewidth=0)

plt.text(xrec - 57., zrec - 22., "Receptor", fontdict=dict(color=crec))
plt.text(xprot - 60., zprot + 14., "Protein", fontdict=dict(color=cprot))

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.axis('off')
plt.title("Solid-state pore")

###############################################################

plt.figure("pugh", figsize=(1.7, 1.7))
pore = get_pore(geoname="pughcyl",
   R=32, Htop=37, Hbot=27)
rw = RandomWalk(pore, rMolecule=2.0779)

rprot = 3.
zprot = rw.ztop + 3.
xprot = -3.

ax = rw.plot_pore()
add_ball(ax, xprot, zprot, rprot, facecolor=cprot, linewidth=0)
plt.text(xprot - 20., zprot + 3., "Protein", fontdict=dict(color=cprot))

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.axis('off')
plt.title("DNA origami pore")

import folders
savefigs("pore", folders.FIGDIR_CWD + "/fig1", ending=".pdf")