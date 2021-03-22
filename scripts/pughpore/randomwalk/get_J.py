import numpy as np
import nanopores.tools.fields as f
import matplotlib.pyplot as plt
import os

HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")


DATADIR = os.path.join(HOME, "fields")
import nanopores.tools.fields as f
f.set_dir(DATADIR)
h = 2.
field = f.get_fields("pughcenter", bulkcon=1e3, Qmol=8, h=h)
z = [x[2] for x in field["x"]]
J = [j*1e12 for j in field["J"]]
I = sorted(list(range(len(z))), key=lambda k: z[k])
z1 = [z[i] for i in I]
J1 = [J[i] for i in I]
from scipy.interpolate import interp1d
Jf = interp1d(z1,J1)
if __name__ == '__main__':
    plt.plot(z1, J1, "s-", label="molecule, h=%s" %h)
    plt.xlabel("z position of molecule [nm]")
    plt.ylabel("current [pA]")
    plt.title("current at -100mV for molecule along pore center")
    plt.legend()
    plt.tight_layout()
    plt.show()
    X=np.linspace(-30.,30.,200)
    plt.plot(X,Jf(X))
    plt.tight_layout()
    plt.show()
