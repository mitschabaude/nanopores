import nanopores.geometries.pughpore as pughpore
import nanopores as nano
#fields.set_dir("/home/benjamin/Dropbox/Paper Howorka/data/fields")

import os

HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")

DATADIR = os.path.join(HOME, "Dropbox", "nanopores", "fields")

import nanopores.tools.fields as f
f.set_dir(DATADIR)

def zsorted(data, field):
    z = [x[2] for x in data["x"]]
    J = data[field]
    I = sorted(range(len(z)), key=lambda k: z[k])
    z1 = [z[i] for i in I]
    J1 = [J[i] for i in I]
    return z1, J1

geop = nano.Params(pughpore.params)
rMolecule = geop.rMolecule
N = 2e4

data = f.get_fields("pugh_diffusivity2D", rMolecule=rMolecule)#, h=4., Nmax=N)
Z, D = zsorted(data, "D")
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
#D=np.load('D.npy')
#Z=np.load('D_Z.npy')
D = np.array(D)
Z = np.array(Z)

a, b = -7.3, 7.3
D_N = np.zeros(Z.shape[0])
for i in range(Z.shape[0]):
    if Z[i] >a and Z[i]<b:
        D_N[i]=np.mean(np.array([D[j] for j in np.arange(i-3,i+4)]))
    else:
        D_N[i]=D[i]

DD = np.array([(D_N[i+1]-D_N[i-1])/(Z[i+1]-Z[i-1]) for i in np.arange(1,Z.shape[0]-1)])
DD = np.concatenate((np.array([0.]),DD))
DD = np.concatenate((DD,np.array([0.])))
f = interp1d(Z,D_N)
fp = interp1d(Z,DD)
if __name__ == "__main__":
    plt.plot(Z,D_N,color='r')
    plt.plot(Z,D,color='b')
    plt.scatter(Z,D)
    plt.scatter(Z,D_N)
    plt.plot(Z,DD,color='green')
    plt.scatter(Z,DD,color='green')
    plt.plot(Z,f(Z),color='k')
    plt.plot(Z,fp(Z),color='k')
    plt.show()
