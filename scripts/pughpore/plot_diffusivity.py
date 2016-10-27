# (c) 2016 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
from diffusion import calculate_diffusivity2D
import nanopores.tools.fields as f
import nanopores

def zsorted(data, field):
    z = [x[2] for x in data["x"]]
    J = data[field]
    I = sorted(range(len(z)), key=lambda k: z[k])
    z1 = [z[i] for i in I]
    J1 = [J[i] for i in I]
    return z1, J1
    
# points
H = 50.
Z = np.linspace(-H, H, 96)
X = [[0.,0.,z] for z in Z]

fig, ax = plt.subplots(figsize=(5, 4))

# get data
for r in [0.152, 0.167, 0.25]:
    data = calculate_diffusivity2D(X, nproc=6, rMolecule=r)
    #data = f.get_fields("pugh_diffusivity2D", rMolecule=0.152, H=100., h=.5)
    Z, D = zsorted(data, "D")

    # plot
    ax.plot(Z, D, "s-", label="r=%.3f, N=270k" %r)
    ax.set_xlabel("z position of molecule [nm]")
    ax.set_ylabel("D/D0")
    ax.set_title("rel. diffusivity (2D model)")
    
# coarser calculation for remaining radii
for r in [0.5, 1., 1.5, 2.0779]:
    N = 2e4
    data = calculate_diffusivity2D(X, nproc=6, rMolecule=r, h=4., Nmax=N)
    #data = f.get_fields("pugh_diffusivity2D", rMolecule=0.152, H=100., h=.5)
    Z, D = zsorted(data, "D")

    # plot
    ax.plot(Z, D, "s-", label="r=%.3f, N=20k" %r)
    ax.set_xlabel("z position of molecule [nm]")
    ax.set_ylabel("D/D0")
    ax.set_title("rel. diffusivity (2D model)")

ax.legend()
from ..howorka.folders import FIGDIR
nanopores.savefigs("pugh_diffusivity", FIGDIR)
plt.show()

