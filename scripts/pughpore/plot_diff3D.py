# (c) 2016 Gregor Mitscha-Baude
import nanopores
from matplotlib import pyplot as plt
from itertools import product
from folders import fields, FIGDIR
#fields.update()

D2D = fields.get_field("pugh_diff2D_test", "D")[0]
data = fields.get_fields("pugh_diff3D_test", bulkbc=True)

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

x = [z[0] for z in data["x"]]
data, x = fields._sorted(data, x)
dstr = ["x", "y", "z"]
for i in range(3):
    Dxx = [D[i][i] for D in data["D"]]
    style = "s-"
    plt.plot(x, Dxx, style, label=r"$D_{%s%s}$" % (dstr[i], dstr[i]))

plt.plot(x, [D2D]*len(x), "--k", label="2D cyl.")
plt.xlabel("distance from pore center [nm]")
plt.ylabel("diffusivity relative to bulk")
plt.legend(loc="best") #bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
plt.gcf().set_size_inches(5, 4)
#nanopores.savefigs("pugh_diff3D", FIGDIR)
plt.show()