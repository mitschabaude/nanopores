# -*- coding: utf-8 -*
# (c) 2018 Gregor Mitscha-Baude
import nanopores.plots as plots

import matplotlib.pyplot as plt
#plt.rc("font", **{"sans-serif": ["cmss10"]})
#plt.rc("font", **{"sans-serif": ["CMU Sans Serif"]})
#from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import numpy as np
import nanopores
import nanopores.geometries.pughpore as pughpore
from nanopores.models.pughpore import polygon as pughpolygon
from nanopores.models.pughpoints import plot_polygon
from nanopores.tools import fields
fields.set_dir_mega()

up = nanopores.Params(pughpore.params)
l0 = up.l0
l1 = up.l1
l2 = up.l2
l3 = up.l3
l4 = up.l4
hpore = up.hpore
hmem = up.hmem
h2 = up.h2
h1 = up.h1
h4 = up.h4

figsize1 = (2.2, 1.7)
figsize2 = (10, 2.21)
lw = .5

# eventspara_nobind_traj_0.00109600_0424_1.7e+07_3.0e+04_1.9e-01_3.0e-0123.0
# eventspara_nobind_traj_0.00236500_0316_1.7e+07_3.0e+04_1.9e-01_3.0e-0123.0
x, y, z, j, t, b1, b2 = fields.get(
        "eventspara_nobind", "X", "Y", "Z", "J", "T", "b1", "b2")

ind = [424, 316] # indices
loc = [8, 9]
colors = [plots.set_sl(plots.colors.pure, 1., .8), plots.colors.mediumdark]

extra_lines = [(0.5*16, -0.5*hpore), (0.5*16, -0.5*hpore + hmem)]

for m in 0, 1:
    i = ind[m]
    X = x[i].load()
    Y = y[i].load()
    Z = z[i].load()
    T = t[i].load()
    J = j[i].load()
    
    I0 = 7.523849e-10
    bind1 = np.where(T>1e6)[0]
    bind2 = np.intersect1d(np.where(T<=1e6)[0], np.where(T>100.)[0])
    amplitude = I0 - np.inner(J, T)/np.sum(T)
    for k in range(1, T.shape[0]):
        T[k] = T[k] + T[k-1]
    tau_off = T[-1]
    #J = J*1e12
    A = (I0 - J)/I0 * 100.
    
    fig=plt.figure(figsize=figsize1)
        
    #color2='#ff0000'
    #color1='#ff9900'
    color3='#aaaaaa'
    
    # #b1 = [[[l1/2.,17.],[l1/2.,19.]],[[l3/2.,-hpore/2.],[l3/2.,hpore/2.-h2],[l2/2.,hpore/2.-h2],[l2/2.,14.]]]
    # for seq in b1:
    #     x= [p[0] for p in seq]
    #     xm=[-p[0] for p in seq]
    #     y= [p[1] for p in seq]
    #     plt.plot(x,y,color=color1,linewidth=2.)
    #     plt.plot(xm,y,color=color1,linewidth=2.)
    # #b2 = [[[l3/2.-.5,-3.],[l3/2.-.5,11.]]]
    # for seq in b2:
    #     x= [p[0] for p in seq]
    #     xm=[-p[0] for p in seq]
    #     y= [p[1] for p in seq]
    #     plt.plot(x,y,color=color2,linewidth=2.)
    #     plt.plot(xm,y,color=color2,linewidth=2.)
    
    plt.plot(X, Z, linewidth=lw, c=colors[m])
    #longer = plt.scatter(X[bind1],Z[bind1],s=50,marker='h',c=color2,linewidth=0.)
    #shorter = plt.scatter(X[bind2],Z[bind2],s=25,marker='h',c=color1,linewidth=0.)
    #start = plt.scatter([X[0]],[Z[0]],s=30,marker='x',c=color3,linewidth=2.)
    
    #patches=[start]
    # labels=['Start']
    
    # if len(bind1)>0:
    #     patches=patches+[longer]
    #     labels+=['Longer bindings']
    # if len(bind2)>0:
    #     patches=patches+[shorter]
    #     labels+=['Shorter bindings']

    c = 1 - 2*m
    plt.annotate("", (X[0], Z[0]),
                xytext=(X[0] + c*4, Z[0] + 4), 
                color=color3, horizontalalignment='left',
                arrowprops=dict(arrowstyle="->", color=color3,
                patchA=None,
                shrinkA=0,
                shrinkB=0))
    plt.text(X[0] + c*4, Z[0] + 5, "Start", color=color3,
        ha=("left" if m == 0 else "right"))
    
    #plt.legend(patches, labels, scatterpoints=1, loc=(.43,.12))
    ax = plt.gca()
    ax.set_aspect("equal")
    
    ax.set_xlim([-20., 60.])
    ax.set_ylim([-25., 35.5])
    
    ax.set_axis_off() # <=> plt.axis("off")
    ax.get_xaxis().set_visible(False) # so that white space disappears!
    ax.get_yaxis().set_visible(False)
    
    #ax.set_xticks([])
    #ax.set_yticks([])
    #plt.axis('off')
    
    plot_polygon(ax, pughpolygon(rmem=100.), linewidth=lw)
    plot_polygon(ax, extra_lines, linewidth=lw)
    
    plt.axes([.6, .45, .3, .3])
    #plt.title('Current  $A/I_0$')
    plt.title("Current signal")
    ax = plt.gca()
    
    if tau_off < 1e3:
        #t = np.linspace(0., tau_off, 3)
        fac = 1.
        #ax.set_xlabel("Time [ns]")
    elif tau_off < 1e6 and tau_off >= 1e3:
        #t = np.linspace(0., tau_off*1e-3, 3)
        fac = 1e-3
        #ax.set_xlabel(u"Time [µs]")
    else:
        #t = np.linspace(0., tau_off*1e-6, 3)
        fac = 1e-6
        #ax.set_xlabel("Time [ms]")
        
    ax.set_xlabel("Time")
        
    T = T*fac
    plt.plot(T, A, color="#000000", linewidth=lw)
    #yt = np.linspace(580., 760, 4)
    plt.ylabel(r"$A/I_0$ [%]", labelpad=1.5)
    ax.set_ylim([-5, 25])
    #plt.yticks([0., 10., 20.], ["0%", "10%", "20%"])
    plt.yticks([0., 10., 20.], [0, 10, 20])
    #plt.xticks([0, 1], ["0", u"1µs"])
    plt.xticks([], [])
    ax.invert_yaxis()
    
    scalebar = AnchoredSizeBar(ax.transData, 1, ur"1$\,$µs", loc[m], 
                   pad=0.2,
                   color="#999999",
                   frameon=False,
                   size_vertical=0.5,
                   fontproperties = fm.FontProperties(size=7)
                   )
    ax.add_artist(scalebar)
    
    #ax.set_xticks(t)
    #xfmt=FormatStrFormatter('%.1f')
    #ax.xaxis.set_major_formatter(xfmt)
    #ax.set_xlim([-4e-2*tau_off*fac,(1.+4e-2)*tau_off*fac])
    
    #fig.tight_layout()
    
##### plot binding types
from matplotlib import transforms

def rainbow_text(x, y, strings, colors, ax=None, **kw):
    """
    Take a list of ``strings`` and ``colors`` and place them next to each
    other, with text strings[i] being shown in colors[i].

    This example shows how to do both vertical and horizontal text, and will
    pass all keyword arguments to plt.text, so you can set the font size,
    family, etc.

    The text will get added to the ``ax`` axes, if provided, otherwise the
    currently active axes will be used.
    """
    if ax is None:
        ax = plt.gca()
    t = ax.transData
    canvas = ax.figure.canvas

    # horizontal version
    for s, c in zip(strings, colors):
        text = ax.text(x, y, s + " ", color=c, transform=t, **kw)
        text.draw(canvas.get_renderer())
        ex = text.get_window_extent()
        t = transforms.offset_copy(text._transform, x=ex.width, units='dots')

def configure_ax(ax):
    plot_polygon(ax, pughpolygon(rmem=100.), linewidth=lw)
    plot_polygon(ax, extra_lines, linewidth=lw)

    ax.set_aspect("equal")
    ax.set_xlim([-20., 20.])
    ax.set_ylim([-33., 40])
    ax.set_axis_off() # <=> plt.axis("off")
    ax.get_xaxis().set_visible(False) # so that white space disappears!
    ax.get_yaxis().set_visible(False)
    
def plot_seq(ax, seq, color, lw):
    x = [p[0] for p in seq]
    xm = [-p[0] for p in seq]
    y = [p[1] for p in seq]
    ax.plot(x, y, color=color, linewidth=lw)
    ax.plot(xm, y, color=color, linewidth=lw)

bshort = pughpolygon(rmem=10)[:7]
blong = fields.get("eventsnew_onlyone_2_", "b2")[0]
blong1 = [(p[0] - .2, p[1]) for p in blong]
blong2 = fields.get("eventsnew_both_1_", "b2")[0]
blong2 = [(p[0] - .5, p[1]) for p in blong2]
    
texty = 32
line = 6
clong = plots.set_sl(plots.colors.experiment, 1., 0.3)
cshort = plots.set_sl(plots.colors.experiment, .7, 0.8)
lw2 = 3*lw

# no binding
fig, ax = plt.subplots(num="bindzones_0", figsize=figsize2)
configure_ax(ax)
plt.text(-11, texty, "No binding", color="k")

# one binding
fig, ax = plt.subplots(num="bindzones_1", figsize=figsize2)
configure_ax(ax)
plt.text(-12.5, texty, "Long", color=clong)
plt.text(-1, texty, "binding", color='k')
plot_seq(ax, blong1, clong, 5*lw)

# two bindings
fig, ax = plt.subplots(num="bindzones_2", figsize=figsize2)
configure_ax(ax)
plt.text(-14.5, texty + .5*line, 'Long', color=clong)
plt.text(-3.3, texty + .5*line, 'and', color='k')
plt.text(5.5, texty + .5*line, 'short', color=cshort)
plt.text(-7.5, texty - .5*line, 'binding', color='k')
plot_seq(ax, bshort, cshort, lw2)
plot_seq(ax, blong2, clong, lw2)


import folders
figdir = folders.FIGDIR_CWD + "/traces"
nanopores.savefigs("trace", figdir, ending=".pdf", pad_inches=0)