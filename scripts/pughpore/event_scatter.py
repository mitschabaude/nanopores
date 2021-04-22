# (c) 2017 Gregor Mitscha-Baude
import nanopores.plots as plots
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
from nanopores.tools import fields, savefigs

fields.set_dir_mega()
t_exp, A_exp = fields.get("events_pugh_experiment", "t", "drop")
t_sim, A_sim = fields.get("eventsnew_both_1_", "t", "a")

plt.figure("pugh", figsize=(1.7, 1.6))
plt.scatter(t_sim, A_sim, s=3, alpha=0.5, linewidth=0,
    color=plots.set_sl(plots.colors.pure, .8, .6), zorder=100, label="Sim.")
plt.scatter(t_exp, A_exp, s=3, alpha=0.5, linewidth=0,
    color=plots.colors.experiment, label="Exp.")
plt.xscale("log")
plt.ylim(-2, 40)
plt.xlim(2e-3, .3e3)

ax = plt.gca()
ax.invert_yaxis()
# ax.set_yticklabels([])
ax.set_yticks([0, 15 ,30])
ax.set_xticks([0.01, 1, 100])
ax.set_xticklabels([0.01, 1, 100])
ax.tick_params(axis="y", direction='in', length=2, pad=-5)
ax.set_yticklabels([0, 15, 30], fontdict={"ha": "left"})
# length=6, width=2, colors='r', grid_color='r', grid_alpha=0.5)

#plt.xlabel(r"$\tau$ off [ms]")
#plt.ylabel(r"$A/I_0$ [%]")
plt.xlabel("Dwell time [ms]")
plt.ylabel("Amplitude [%]")

plots.removeTopRightFrame()

plt.legend(frameon=True, scatterpoints=3,
    loc="upper right", bbox_to_anchor=(1.06, 1.0)) #, loc=(0.01, 0.01))
plt.title("Current events")

import folders
savefigs("events", folders.FIGDIR_CWD + "/fig1", ending=".pdf")