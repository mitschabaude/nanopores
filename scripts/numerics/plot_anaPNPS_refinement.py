import nanopores
import matplotlib.pyplot as plt
from mysolve import load_estimators as load

figsize = 5*0.8, 4*0.8

est2D = load("anaPNPS_2D")
est3D = load("anaPNPS_3D")

fig = plt.figure("refinement", figsize=figsize)
est2D["Jvol"].name = r"$|J_h - J|/J$ (2D)"
est3D["Jvol"].name = r"$|J_h - J|/J$ (3D)"
est2D["Jvol"].plot(rate=-1., fig=False)
est3D["Jvol"].plot(rate=-2./3, fig=False)

ax = plt.gca()
ax.lines[1].set_label("")
ax.lines[3].set_label("")

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

from folders import FIGDIR
nanopores.savefigs("anaPNPS", FIGDIR)