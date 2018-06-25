# (c) 2017 Gregor Mitscha-Baude
import matplotlib.pyplot as plt
from nanopores.tools import fields, savefigs

fields.set_dir_mega()
t_exp, A_exp = fields.get("events_pugh_experiment", "t", "drop")
t_sim, A_sim = fields.get("eventsnew_both_1_", "t", "a")

plt.figure("pugh", figsize=(2.7, 2.3))
plt.scatter(t_sim, A_sim, 4, color="C0", zorder=100, label="Simulation")
plt.scatter(t_exp, A_exp, 4, color="r", label="Experiment")

plt.ylim(ymax=51)
plt.xlim(.5e-5, 3e2)
plt.xscale("log")
plt.xlabel(r"$\tau$ off [ms]")
plt.ylabel(r"$A/I_0$ [%]")
plt.gca().invert_yaxis()
plt.legend(frameon=False, scatterpoints=3, loc=(0.01, 0.01))
plt.title("Current events")

import folders
savefigs("output", folders.FIGDIR + "/fig1", ending=".eps")