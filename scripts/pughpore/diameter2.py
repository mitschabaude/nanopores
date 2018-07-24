# -*- coding: utf-8 -*
# (c) 2017 Gregor Mitscha-Baude
from matplotlib import rcParams, rc
rcParams.update({
    "font.size" : 7,
    "axes.titlesize" : 7,
    "font.family" : "sans-serif",
    "font.sans-serif" : ["Helvetica"],
    "lines.linewidth" : 1,
    "lines.markersize" : 5,
})
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.collections as mcol
from matplotlib.legend_handler import HandlerTuple, HandlerLineCollection
from matplotlib.lines import Line2D
#import matplotlib.ticker as ticker
import nanopores
import nanopores.models.pughpore as pugh
from nanopores.models.diffusion_interpolation import (cache_pugh_diffusivity,
                                                      diff_profile_z_pugh)
from nanopores.tools import fields
fields.set_dir_mega()

#### color code
color_open = "C0" #"#0066ff"
color_open_light = mcolor.to_rgb(color_open) + (0.3,)
color_closed = "#00cc00"
color_closed_light = mcolor.to_rgb(color_closed) + (0.3,)
color_drop = "red"
color_drop_light = mcolor.to_rgb(color_drop) + (0.3,)

msize = 4.
msize2 = 4.2

figsize1 = (1.73, 1.53)
figsize2 = (1.73, 1.37)

lc = mcol.LineCollection(2 * [[(0, 0)]], linestyles=["--", "--"],
                        colors=[color_open, color_closed])

# mpl hack for legend
class HandlerDashedLines(HandlerLineCollection):
    "Custom Handler for LineCollection instances."
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        # figure out how many lines there are
        numlines = len(orig_handle.get_segments())
        xdata, xdata_marker = self.get_xdata(legend, xdescent, ydescent,
                                             width, height, fontsize)
        leglines = []
        # divide the vertical space where the lines will go
        # into equal parts based on the number of lines
        ydata = ((height) / (numlines + 1)) * np.ones(xdata.shape, float)
        # for each line, create the line at the proper location
        # and set the dash pattern
        for i in range(numlines):
            legline = Line2D(xdata, ydata * (numlines - i) - ydescent)
            self.update_prop(legline, orig_handle, legend)
            # set color, dash pattern, and linewidth to that
            # of the lines in linecollection
            try:
                color = orig_handle.get_colors()[i]
            except IndexError:
                color = orig_handle.get_colors()[0]
            try:
                dashes = orig_handle.get_dashes()[i]
            except IndexError:
                dashes = orig_handle.get_dashes()[0]
            try:
                lw = orig_handle.get_linewidths()[i]
            except IndexError:
                lw = orig_handle.get_linewidths()[0]
            if dashes[0] is not None:
                legline.set_dashes(dashes[1])
            legline.set_color(color)
            legline.set_transform(trans)
            legline.set_linewidth(lw)
            leglines.append(legline)
        return leglines

dparams = {2: dict(Nmax=1e5, dim=2, rMolecule=0.11, h=1.0),
           3: dict(Nmax=2e6, dim=3, rMolecule=0.11, h=2.0)}

# I for different pore diameters
@pugh.solvers.cache_forcefield("pugh_Idiam4_", pugh.defaultp)
def Idiam(diam, **params):
    dim = params["dim"]
    x0 = params.pop("x0")
    params.pop("diamPore")
    Jon = []
    Joff = []

    for d in diam:
        # get diffusivity interpolation
        cache_pugh_diffusivity(geoname="pugh2", diamPore=d, **dparams[dim])

        # current WITHOUT molecule
        setup = pugh.Setup(x0=None, diamPore=d, **params)
        pb, pnps = pugh.solve(setup, visualize=True)
        Jon.append(pnps.evaluate(setup.phys.CurrentPNPS)["J"])

        # current WITH molecule
        setup = pugh.Setup(x0=x0, diamPore=d, **params)
        pb, pnps = pugh.solve(setup, visualize=True)
        Joff.append(pnps.evaluate(setup.phys.CurrentPNPS)["J"])

    return dict(Jon=Jon, Joff=Joff)

params = {2: dict(dim=2, h=1., Nmax=1e5, x0=[0.,0.,0.], bV=-0.08,
                  diffusivity="Dpugh2"),
          3: dict(dim=3, h=2., Nmax=6e5, x0=[0.,0.,0.], bV=-0.08,
                  diffusivity="Dpugh2", stokesiter=True, cheapest=False)}

diam = {2: [4.18, 4.25, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.93, 5.2, 5.5, 6., 6.5, 7., 7.5],
        3: [4.22, #4.3, 
            4.4, #4.5, 
            4.6, 4.8, 5., 5.5, 6., 6.5, 7., 7.5]}
        #3: [4.18, 4.2, 4.25, 4.28, 4.3, 4.35, 4.4, 4.5, 4.555, 4.6, 4.8, 5., 5.5, 6., 6.5, 7., 7.5, 8.]}

fit = {2: 4.93, 3: 4.4}

calc = nanopores.user_param(calc=False)

if calc:
    for dim in 2, 3:
        for d in diam[dim]:
            diff_profile_z_pugh(nproc=6, diamPore=d)
#diam = {2: [3.7, 3.8, 3.9, 4.0, 4.141, 4.4, 4.6,
#            4.8, 5., 5.5, 6., 6.65],# 7.], # 7.5, 8.],
#        3: [4.18, 4.23, 4.3, 4.4,
#            4.5, 4.6, 4.8, 5.2, 5.5, 6., 7., 7.501]}

# semi-bad values: 4.17, 4.19, 4.21, 4.22, 4.25, 4.275, 4.7,
# bad values: 4.19, 4.9, 4.995, 5.095,

# for old setup:
#        3: [4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 4.23, 4.25, 4.275, 4.3, 4.4,
#            4.45, 4.5, 4.6, 4.8, 5., 5.5, 6., 6.5, 7., 7.5, 8.]}

for dim in 2, 3:
    fig = plt.figure("abs_%dD" % dim, figsize=figsize1)
    result = Idiam(diam[dim], calc=calc, nproc=2, **params[dim])
    d = result["x"]
    print "diameters (%dD):" % dim, d
    print "missing:", set(diam[dim]) - set(d)
    
    xlabel = {2: "Pore diameter [nm]", 3: "Channel width [nm]"}[dim]

    n = len(d)
    Jon = 1e9*np.array(result["Jon"])
    Joff = 1e9*np.array(result["Joff"])
    sim1 = plt.plot(d, Jon, "s", markersize=msize, color=color_open_light,
             mec="None", mfc=color_open, label="No protein")
    plt.plot(d, Jon, "-", color=color_open_light)
    exp1 = plt.plot([0,10], [2.29*0.08]*2, "--", color=color_open, label="Exp.")
    plt.fill_between([0,10], [(2.29 - 0.26)*0.08]*2,
                     [(2.29 + 0.26)*0.08]*2,
                     color=color_open, alpha=0.3, lw=0)
                     
    sim2 = plt.plot(d, Joff, "v", markersize=msize2, color=color_closed_light,
             mec="None", mfc=color_closed, label="With protein")
    plt.plot(d, Joff, "-", color=color_closed_light)
    exp2 = plt.plot([0,10], [2.29*0.08*(1 - 0.262)]*2, "--", 
                     color=color_closed, label="Exp.")
    plt.fill_between([0,10], [(2.29*(1 - 0.262) - 0.26)*0.08]*2,
                        [(2.29*(1 - 0.262) + 0.26)*0.08]*2,
                        color=color_closed, alpha=0.3, lw=0)
    plt.xlabel(xlabel)
    plt.ylabel("Current [nA]")

    plt.xlim(4., 7.6)
    plt.ylim(0., 1.250)
    #locx, _ = plt.xticks()
    plt.yticks([0, 0.5, 1], [0, 0.5, 1])
    plt.gca().set_yticks([.1,.2,.3,.4,.6,.7,.8,.9,1.1,1.2], minor=True)
    plt.axvline(x=2*2.0779, linestyle="--", color="#666666")
    plt.annotate(u"Protein ⌀", (2*2.0779, 1.15),
                 xytext=(2*2.0779 + 0.6, 1.15-.02), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))
    plt.text(4.5, .82, "W/o protein", color=color_open)
    plt.text(6., .58, "W/ protein", color=color_closed)
    #plt.title("influence of pore diameter on current (%dD)" %dim)
    plt.legend([(sim1[0], sim2[0]), lc],
                ["Sim.", "Exp."],
    #plt.legend([sim1[0], sim2[0], lc],
               #["Sim. w/o prot.", "Sim. w/ prot.", "Experiment"],
               #["No prot.", "With prot.", "Exp."],
                handler_map={tuple: HandlerTuple(ndivide=None),
                             type(lc): HandlerDashedLines()},
                loc="lower right",
                #borderaxespad=0.1, borderpad=0.3,
                #handletextpad=0.4
                )
    #plt.legend(loc="lower right")
    #fig.tight_layout()

    fig = plt.figure("drop_%dD" % dim, figsize=figsize2)
    drop = (1. - Joff/Jon)*100
    plt.plot(d, drop, "o", markersize=msize, color=color_drop, label="Sim.")
    plt.plot(d, drop, "-", color=color_drop_light)
    plt.plot([0,10], [26.2]*2, "--r", label="Exp.")
    plt.fill_between([0,10], [26.2 - 0.7]*2, [26.2 + 0.7]*2,
                     color=color_drop, alpha=0.3, lw=0)
    plt.xlim(4., 7.6)
    #plt.xticks([4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5])
    #plt.ylim(ymin=0., ymax=50.)
    #plt.xlabel("pore diameter [nm]")
    plt.xlabel(xlabel)
    plt.ylabel(r"$A/I_0$ [%]")
#    else:
#        loc, _ = plt.yticks()
#        plt.yticks(loc, [])
    ymax = 100. if dim==2 else 40.
    htext = ymax*0.9
    plt.ylim(ymin=0., ymax=ymax)

    #loc, _ = plt.xticks()
    #plt.xticks(locx, [])
    #plt.yscale("log")
    #ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*0.01))
    #plt.gca().yaxis.set_major_formatter(ticks)

    plt.gca().invert_yaxis()
    plt.axvline(x=fit[dim], ymin=0.,
                ymax=(ymax - 26.2)/ymax,
                color=color_drop_light, zorder=-90)
    plt.scatter([fit[dim]], [26.2], s=200, c=color_drop_light, linewidths=0)
    plt.annotate("%.1f nm" % (fit[dim],), (fit[dim], 26.2),
             xytext=(fit[dim] + 0.16, 26.2 + 5.5*(ymax/40.)),
             color=color_drop_light)

    plt.axvline(x=2*2.0779, linestyle="--", color="#666666")
    #htext = 45 # 2
    xofftext = 0.5 # 0.2
    
#    plt.annotate(u"Protein ⌀", (2*2.0779, htext),
#             xytext=(2*2.0779 + xofftext, htext-1.), color="#666666",
#             arrowprops=dict(arrowstyle="->", color="#666666"))

#    plt.annotate("diameter of trypsin", (2*2.0779, 2),
#                 xytext=(2*2.0779 + 0.2, 2-0.2), color="#666666",
#                 arrowprops=dict(arrowstyle="->", color="#666666"))
#
    if dim==3:
        plt.legend(loc="center right", bbox_to_anchor=(1., .58))
    else:
        plt.legend(loc="lower right")
    #fig.tight_layout()
    #plt.subplots_adjust(left=0.2)
    #plt.legend(loc="center right" if dim==3 else "best")
    #fig.savefig(DIR + name + "_" + label + ".eps", bbox_inches=bbox_inches)


import folders
pugh.nano.savefigs("new", folders.FIGDIR + "/Idiam2", ending=".pdf")
#plt.show()