# (c) 2016 Gregor Mitscha-Baude
import nanopores
from matplotlib import pyplot as plt
from itertools import product
from folders import fields, FIGDIR
import numpy as np
#fields.update()

#D2D = fields.get_field("pugh_diff2D_test", "D")[0]
data = fields.get_fields("pugh_diff3D_cross", bulkbc=True, rMolecule=2.0779)

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
def smooth3(l):
    A=np.array(l)
    B=A[:]
    ker=np.array([1./3,1./3,1./3])
    n=int(ker.shape[0]/2.)
    for i in range(n,A.shape[0]-n):
        B[i]=np.inner(A[i-n:i+n+1],ker)
    return list(B)
def smooth5(l):
    A=np.array(l)
    B=A[:]
    ker=np.array([.2,.2,.2,.2,.2])
    n=int(ker.shape[0]/2.)
    for i in range(n,A.shape[0]-n):
        B[i]=np.inner(A[i-n:i+n+1],ker)
    return list(B)





x = [z[0] for z in data["x"]]
data, x = fields._sorted(data, x)
dstr = ["x", "y", "z"]
x.append(1.0)

Dxx = [D[0][0] for D in data["D"]]
Dyy = [D[1][1] for D in data["D"]]
Dzz = [D[2][2] for D in data["D"]]
Dxx.append(0.)
Dyy.append(0.)
Dzz.append(0.)
style = "s-"
plt.plot(x, smooth5(Dxx), style, label=r"$D_{%s%s}$" % (dstr[0], dstr[0]))
plt.plot(x, smooth3(Dyy), style, label=r"$D_{%s%s}$" % (dstr[1], dstr[1]))
plt.plot(x, smooth3(Dzz), style, label=r"$D_{%s%s}$" % (dstr[2], dstr[2]))

#plt.plot(x, [D2D]*len(x), "--k", label="2D cyl.")
plt.xlabel("distance from pore center [nm]")
plt.ylabel("diffusivity relative to bulk")
plt.legend(loc="lower left") #bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
plt.gcf().set_size_inches(5, 4)
#nanopores.savefigs("pugh_diff3D_r0.11", FIGDIR)
plt.show()
