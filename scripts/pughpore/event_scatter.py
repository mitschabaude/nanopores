# (c) 2017 Gregor Mitscha-Baude
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
rcParams.update({
    "font.size" : 7,
    "axes.titlesize" : 7,
    "font.family" : "sans-serif",
    "font.sans-serif" : ["Helvetica"],
    "lines.linewidth" : 1,
    "lines.markersize" : 5,
})
from nanopores.tools import fields, savefigs

fields.set_dir_mega()
t_exp, A_exp = fields.get("events_pugh_experiment", "t", "drop")
t_sim, A_sim = fields.get("eventsnew_both_1_", "t", "a")

plt.figure("pugh", figsize=(1.7, 1.6))
plt.scatter(t_sim, A_sim, s=3, alpha=0.5, linewidth=0, color="C0", zorder=100, label="Simulation")
plt.scatter(t_exp, A_exp, s=3, alpha=0.5, linewidth=0, color="r", label="Experiment")
plt.xscale("log")
plt.ylim(ymax=51)
plt.xlim(.5e-5, 3e2)

ax = plt.gca()
ax.invert_yaxis()
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_yticks([])
ax.set_xticks([])

#plt.xlabel(r"$\tau$ off [ms]")
#plt.ylabel(r"$A/I_0$ [%]")
plt.xlabel("Dwell time")
plt.ylabel("Amplitude")

plt.legend(frameon=False, scatterpoints=3, loc=(0.01, 0.01))
plt.title("Current events")

import folders
savefigs("output", folders.FIGDIR + "/fig1", ending=".pdf")