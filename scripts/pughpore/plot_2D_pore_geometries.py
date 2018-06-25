# (c) 2018 Gregor Mitscha-Baude
"plot simulated pores schematically with protein/receptors"
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from nanopores import get_pore, savefigs
from nanopores.models.randomwalk import RandomWalk, Ball

# for plotting proteins andreceptors
ball_settings = dict(facecolor="#eeeeee", edgecolor="k", linewidth=1.)
def add_ball(ax, x, y, r, **settings):
    p = mpatches.Circle((x, y), r, **dict(ball_settings, **settings))
    p.set_zorder(200)
    ax.add_patch(p)

###############################################################

plt.figure("ahem", figsize=(3, 3))
pore = get_pore(geoname="alphahem", R=8, Htop=3.5, Hbot=12.5)
rw = RandomWalk(pore)
ax = rw.plot_pore()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.title("Protein pore")

###############################################################

plt.figure("wei", figsize=(3, 3))
pore = get_pore(geoname="wei",
    R=100, Htop=120, Hbot=80)
rw = RandomWalk(pore)

rrec = 5
zrec = rw.zbot + rrec + (rw.ztop - rw.zbot - 2.*rrec)*0.6
xrec = pore.radius_at(zrec) - rrec - 1.
crec = "red"

rprot = 10
zprot = rw.ztop + 30.
xprot = -20.
cprot = "blue"

ax = rw.plot_pore()
add_ball(ax, xrec, zrec, rrec, facecolor=crec, linewidth=0)
add_ball(ax, xprot, zprot, rprot, facecolor=cprot, linewidth=0)

plt.text(xrec - 45., zrec - 18., "Receptor", fontdict=dict(color=crec))
plt.text(xprot - 55., zprot + 10., "Protein", fontdict=dict(color=cprot))

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.title("Solid-state pore")

###############################################################

plt.figure("pugh", figsize=(3, 3))
pore = get_pore(geoname="pughcyl",
   R=32, Htop=37, Hbot=27)
rw = RandomWalk(pore)

rprot = 3.
zprot = rw.ztop + 5.
xprot = -5.
cprot = "blue"

ax = rw.plot_pore()
add_ball(ax, xprot, zprot, rprot, facecolor=cprot, linewidth=0)
plt.text(xprot - 18., zprot + 1., "Protein", fontdict=dict(color=cprot))

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.title("DNA origami pore")

import folders
savefigs("pore", folders.FIGDIR + "/fig1", ending=".eps")