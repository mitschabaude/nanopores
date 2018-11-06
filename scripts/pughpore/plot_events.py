# -*- coding: utf-8 -*
# (c) 2018 Gregor Mitscha-Baude
import numpy as np
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
from matplotlib import patches as mpatch
from matplotlib.ticker import FormatStrFormatter
from matplotlib.transforms import ScaledTranslation

import nanopores
from nanopores import fields
figdir = nanopores.dirnames.DROPBOX_FIGDIR
fields.set_dir_mega()
a_exp, t_exp = fields.get("events_pugh_experiment", "drop", "t")

# colors and settings
figsize = (2.6, 1.58)

long_bind_color = "C0"

settings = dict(s=3., alpha=0.5, linewidth=0)
s_experiment = dict(settings, facecolors="#888888", label="Experiment", zorder=0)
s_success = dict(settings, facecolors="green", label="Translocation", zorder=2)
s_fail = dict(settings, facecolors="red", label="No transloc.", zorder=1)

def format_ticks(ax):
    ax.set_xscale("log")
    #ax.set_xlabel(u'τ off [µs]')
    ax.set_xlabel(ur"$\tau_\mathrm{off}$ [ms]")
    ax.set_ylabel(r"$A/I_0$ [%]")
    
    ax.set_xlim([5e-5, 2e2])
    ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2])
    #ax.set_xticklabels([ur"100ns", u"1µs", u"10µs", u"100µs",
    #                    u"1ms", u"10ms", "100ms"])
    #xfmt = FormatStrFormatter("%g")
    #ax.xaxis.set_major_formatter(xfmt)
    
    ax.set_ylim([0, 40])
    ax.set_yticks([0., 10., 20., 30., 40.])
    ax.invert_yaxis()

########## simulations without binding ##########
plt.figure("0bind", figsize=figsize)
a, t, ood = fields.get("eventspara_nobind", "a", "t", "ood") 
a, t, success = np.array(a), np.array(t), np.array(ood) == 0

plt.scatter(t[~success], a[~success], **s_fail)
plt.scatter(t[success], a[success], **s_success)
plt.scatter(t_exp, a_exp, **s_experiment)
ax = plt.gca()
format_ticks(ax)

plt.legend(loc="lower left", scatterpoints=3, frameon=False,
           borderaxespad=0.2, handletextpad=0.4)

# encircle long bindings
offset = ScaledTranslation(18, 27, ax.transScale)
tform = offset + ax.transLimits + ax.transAxes
ell = mpatch.Ellipse((0, 0), 1.8, 10, angle=0.0, transform=tform,
                     fc="None", ec=long_bind_color)
ax.add_patch(ell)
plt.text(15, 21, "Assumed\n long bindings", color=long_bind_color,
         horizontalalignment="center")

########## simulations with one binding ##########
plt.figure("1bind", figsize=figsize)
a, t, ood = fields.get("eventsnew_onlyone_2_", "a", "t", "ood") 
a, t, success = np.array(a), np.array(t), np.array(ood) == 0

plt.scatter(t[~success], a[~success], **s_fail)
plt.scatter(t[success], a[success], **s_success)
plt.scatter(t_exp, a_exp, **s_experiment)
format_ticks(plt.gca())

########## simulations with two bindings ##########
plt.figure("2bind", figsize=figsize)
a, t, ood = fields.get("eventsnew_both_1_", "a", "t", "ood") 
a, t, success = np.array(a), np.array(t), np.array(ood) == 0

plt.scatter(t[~success], a[~success], **s_fail)
plt.scatter(t[success], a[success], **s_success)
plt.scatter(t_exp, a_exp, **s_experiment)
format_ticks(plt.gca())

figdir = nanopores.dirnames.HOME + "/Dropbox/Paper Howorka/figures/"
nanopores.savefigs("rw/events", figdir, ending=".pdf")