# (c) 2016 Gregor Mitscha-Baude
import nanopores
from matplotlib import pyplot as plt
from itertools import product
from folders import fields, FIGDIR
#fields.update()
r = 0.11
#D2D = fields.get_field("pugh_diff2D_test", "D")[0]
#data = fields.get_fields("pugh_diff3D_cross", bulkbc=True, rMolecule=2.0779)
data = fields.get_fields("pugh_diff3D", rMolecule=r)

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
fields.set_dir_default()
X, D = fields.get("diffz_pugh", "x", "D", diamPore=6.)
zmin = min([x1[2] for x1 in X], key=lambda x: abs(x))
i = X.index([0., 0., zmin])
D0 = D[i]
Dxx1 = Dxx[-1]
Dzz1 = Dzz[-1]
xlin = linspace(r+1e-3, 3., 100)
dn = [Dn_plane(t, r, N=20) for t in xlin]
dn = [d*D0/dn[-1] for d in dn]
plt.plot(xlin, dn, "-b")
dt = [Dt_plane(t, r) for t in xlin]
dt = [d*D0/dt[-1] for d in dt]
plt.plot(xlin, dt, "-g")
plt.xlim(0., 2.5)

plt.plot(x, Dxx, "ob", label=r"$D_{xx}$")
plt.plot(x, Dyy, "sg", label=r"$D_{yy}$")
plt.plot(x, Dzz, ".r", label=r"$D_{zz}$")




#plt.plot(x, [D2D]*len(x), "--k", label="2D cyl.")
plt.xlabel("x distance from pore wall [nm]")
plt.ylabel("diffusivity relative to bulk")
plt.ylim(0, 1)
plt.axvline(x=0.11, linestyle="--", color="#666666")
plt.annotate("ion radius", (0.11, 0.94),
                 xytext=(0.25, 0.94-0.002), color="#666666",
                 arrowprops=dict(arrowstyle="->", color="#666666"))
plt.yticks([i/10. for i in range(11)])
plt.legend(loc="lower right") #bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
plt.gcf().set_size_inches(4.5, 4.5)
nanopores.savefigs("pugh_Dions", FIGDIR)
plt.show()
