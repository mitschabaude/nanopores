import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerPatch
from collections import OrderedDict
import numpy as np
import nanopores
import os

def plot_circle(ax, x, y, color):
    area = 380    
    circle = ax.scatter(x, y, s=area, c=color)
    #edgewidth = 0
    #edgecolor = "k" # black      
    #circle = ax.plot(x, y, "o", c=color, markersize=10)
    #circle= mpatches.Circle((x, y), radius, fc=color, ec=edgecolor, lw=edgewidth)
    #ax.add_patch(circle)
    return circle


def make_legend_ellipse(legend, orig_handle, xdescent, ydescent,
    width, height, fontsize):
    return mpatches.Ellipse(xy=(0.5*width-0.5*xdescent, 0.5*height-0.5*ydescent),
        width = width+xdescent, height=(height+ydescent))


data = nanopores.load_stuff("robustness_fixedpoint")
data = [np.array(d) for d in data]

bVs = [.01, .02, .05, .1, .2, .5, 1., 2.]
qs = [0.1, 0.25, 0.5, 1., 2.]

X, Y = bVs, qs
xlim = [0.005, 4.]
ylim = [0.05, 4.]

DIR = os.path.expanduser("~") + "/papers/pnps-numerics/figures/"

# --- plot engine ---

def plot_robustness(ax, fpdata, fplabels, fpcolors):    
    handles = OrderedDict() # for legend
    
    for data, label, color in zip(fpdata, fplabels, fpcolors):
        handle = []
        
        # nice painless way to iterate array:
        array = np.nditer(data, flags=['multi_index'])    
        for boolean in array:
            if boolean:
                i, j = array.multi_index
                x, y = X[i], Y[j]
                circle = plot_circle(ax, x, y, color)
                handle.append(circle)
                
        handles[label] = tuple(handle)
    
    # TODO: axis limits
    #ax.set_aspect('equal', 'datalim')
    #ax.autoscale_view(True,True,True)
    ax.set_xscale("log")#, subsx=[2,5])
    ax.set_yscale("log")#, subsy=[2,5])
    
    myhandles = [plt.Rectangle((0, 0), 1, 1, fc=col) for col in fpcolors]
    #ax.legend([h[0] for h in handles.values()], handles.keys(),
    ax.legend(myhandles, list(handles.keys()),
        #bbox_to_anchor=(0.5, 1.05), loc="lower center", borderaxespad=0.,)
        bbox_to_anchor=(0., 1.5), loc="upper left", borderaxespad=0.,)
        #handler_map={mpatches.Circle:HandlerPatch(patch_func=make_legend_ellipse)})
    ax.set_xlabel("voltage bias [V]")
    ax.set_ylabel(r"surface charge density [q/nm$^2$]")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

# --- fixed point ---
fpdata = [
    data[0],
    np.logical_and(data[1], np.logical_not(data[0])),
    np.logical_and(data[2], np.logical_not(data[1])),
    np.logical_not(data[2])
    ]
fplabels = [
    "fixed point",
    "+ voltage schedule (20 iterations)",
    "+ voltage schedule (100 iterations)",
    "no convergence"
    ]
fpcolors = ["#0000ff", "#6666ff", "#ccccff", "w"]

# --- hybrid ---
hydata = [
    data[3],
    np.logical_and(data[4], np.logical_not(data[3])),
    np.logical_not(data[4])
    ]
hylabels = [
    "hybrid",
    "hybrid + PB initial guess",
    "no convergence"
    ]
hycolors = ["#ff0000", "#ffaaaa", "w"]

# --- plot ---

ax = plt.subplot(111)
plot_robustness(ax, fpdata, fplabels, fpcolors)
fig = plt.gcf()
fig.set_size_inches(5.2, 3.5)
plt.savefig(DIR + "robust_fp.eps", bbox_inches='tight')

plt.figure()
ax = plt.subplot(111)
plot_robustness(ax, hydata, hylabels, hycolors)
fig = plt.gcf()
fig.set_size_inches(5.2, 3.5)
plt.savefig(DIR + "robust_hy.eps", bbox_inches='tight')

#plt.show()