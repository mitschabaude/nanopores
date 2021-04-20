# (c) 2016 Gregor Mitscha-Baude
import nanopores
import nanopores.plots as plots
import numpy as np
from matplotlib import pyplot as plt
from itertools import product
from folders import fields, FIGDIR, FIGDIR_CWD

import matplotlib.collections as mcol
from matplotlib.legend_handler import HandlerTuple, HandlerLineCollection
from matplotlib.lines import Line2D

#fields.update()
r = 0.11
#D2D = fields.get_field("pugh_diff2D_test", "D")[0]
#data = fields.get_fields("pugh_diff3D_cross", bulkbc=True, rMolecule=2.0779)
data = fields.get_fields("pugh_diff3D_backup", rMolecule=r)

#def _sorted(data, key):
#    I = sorted(range(len(key)), key=lambda k: key[k])
#    return {k: [data[k][i] for i in I] for k in data}, [key[i] for i in I]

#x = [z[0] for z in data["x"]]
#print len(x)
#data, x = _sorted(data, x)
#dstr = ["x", "y", "z"]
#for i, j in product(range(3), range(3)):
#    Dxx = [D[i][j] for D in data["D"]]
#    style = "s-" if i==j else "--"
#    plt.plot(x, Dxx, style, label=r"$D_{%s%s}$" % (dstr[i], dstr[j]))
#
#plt.plot(x, [D2D]*len(x), "-k", label="2D ref.")
#plt.legend(loc="best")
##plt.legend(bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
##nanopores.savefigs("pugh_diff3D", folders.FIGDIR)
#plt.show()


x = [z[0]-r for z in data["x"]]
data, x = fields._sorted(data, x)
dstr = ["x", "y", "z"]
#print x
x0 = [0.009999999999999662, 0.08555555555555504, 0.16111111111111087,
     0.16157894736842093, 0.23666666666666625, 0.31222222222222207,
     0.31315789473684175, 0.38777777777777744, 0.46333333333333326,
     0.464736842105263, 0.5388888888888886, 0.6144444444444445,
     0.6163157894736843, 0.6899999999999998, 0.7678947368421051,
     0.9194736842105261, 1.071052631578947, 1.2226315789473683,
     1.3742105263157893, 1.5257894736842104, 1.6773684210526314,
     1.8289473684210524, 1.9805263157894732, 2.1321052631578947,
     2.283684210526316, 2.435263157894737, 2.586842105263158,
     2.738421052631579, 2.89]

x = [0.009999999999999662, 0.08555555555555504,
     0.16157894736842093, 0.23666666666666625,
     0.31315789473684175, 0.38777777777777744,
     0.464736842105263, 0.5388888888888886,
     0.6163157894736843, 0.6899999999999998, 0.7678947368421051,
     0.9194736842105261, 1.071052631578947, 1.2226315789473683,
     1.3742105263157893, 1.5257894736842104, 1.6773684210526314,
     1.8289473684210524, 1.9805263157894732, 2.1321052631578947,
     2.283684210526316]

DD = data["D"]
Dxx = [D[0][0] for D in DD if x0[DD.index(D)] in x]
Dyy = [D[1][1] for D in DD if x0[DD.index(D)] in x]
Dzz = [D[2][2] for D in DD if x0[DD.index(D)] in x]

x = [t+r for t in x]

from nanopores.models.diffusion_interpolation import Dn_plane, Dt_plane
from numpy import linspace

fig, ax = plt.subplots(figsize=(1.73, 1.53))
colors = plots.colors
cx = colors.mediumlight
cy = colors.mediumdark
cz = colors.lightpink

Dxx1 = max(Dxx)
Dzz1 = max(Dzz)
xlin = linspace(r+1e-3, 3., 100)
dn = [Dn_plane(t, r, N=20) for t in xlin]
dn = [d*Dxx1/dn[-1] for d in dn]

dt = [Dt_plane(t, r) for t in xlin]
dt = [d*Dzz1/dt[-1] for d in dt]
plt.xlim(0., 1.5)
plt.xticks([0, 0.5, 1., 1.5])

plt.axvline(x=0.11, linestyle="--", color="#666666")

plt.plot(xlin, dn, "-", color=cx, label=r"D$_{xx}$ (r-dep.)")
plt.plot(xlin, dt, "-", color=cy, label=r"D$_{yy}$ (r-dep.)")

dysym = plt.plot(x, Dyy, "o", label=r"D$_{yy}$", color=cy, markersize=3)
dxsym = plt.plot(x, Dxx, "o", label=r"D$_{xx}$", color=cx, markersize=3)
#dzsym = plt.plot(x, Dzz, ".r", label=r"$D_{zz}$", color=cz)


#plt.plot(x, [D2D]*len(x), "--k", label="2D cyl.")
plt.xlabel("Ion distance from pore wall [nm]")
plt.ylabel("Rel. diffusivity")
plt.ylim(0, 1)
plt.xlim(xmax=1.53)
# plt.annotate("Ion radius", (0.11, 0.94),
#                xytext=(0.32, 0.94+0.002), color="#666666",
#                arrowprops=dict(arrowstyle="->", color="#666666"))
plt.yticks([0, 0.5, 1])
plots.removeTopRightFrame(ax)

plt.text(.26, .88, r"D$_{\parallel}$", color=cy, size=8)
plt.text(.8, .66, r"D$_{\bot}$", color=cx, size=8)
#plt.text(.53, .92, r"D$_{zz}$", color=cz)

# legend!!
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

lc = mcol.LineCollection(2 * [[(0, 0)]], linestyles=["-", "-"],
     colors=[cy, cx])

plt.legend([(dysym[0], dxsym[0]), lc], ["LRNH", "r- and z-dep."],
          handler_map={tuple: HandlerTuple(2),
                    type(lc): HandlerDashedLines()},
          loc="lower right",
          frameon=False
          #borderaxespad=0.1, borderpad=0.3,
          #handletextpad=0.4
)

#plt.legend(loc="lower right", frameon=False) #bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)

nanopores.savefigs("pugh/Dions", FIGDIR_CWD, ending=".pdf")
plt.show()
