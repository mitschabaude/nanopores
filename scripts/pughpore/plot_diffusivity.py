# (c) 2016 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
from diffusion import calculate_diffusivity2D
import nanopores.tools.fields as f
import nanopores
import folders

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

fig_big, ax_big = plt.subplots(figsize=(10, 8), num="all")
fig_small, ax_small = plt.subplots(figsize=(6, 4), num="small")
# get data
for r in [0.152, 0.167, 0.25, 2.0779]:
    
    #data = calculate_diffusivity2D(X, nproc=6, rMolecule=r, h=.6, Nmax=2.7e5)
    #data = calculate_diffusivity2D(X, nproc=6, rMolecule=r)
    data = f.get_fields("pugh_diffusivity2D", rMolecule=r, h=.6, Nmax=2.7e5)
    Z, D = zsorted(data, "D")

    # plot
    ax = ax_big
    ax.plot(Z, D, ".-", label="r=%.3f, N=270k" %r)
    ax.set_xlabel("z position of molecule [nm]")
    ax.set_ylabel("D/D0")
    ax.set_title("rel. diffusivity (2D model)")
    
    if r>=0.25: continue
    names = {0.152: r"$\rm{Na}^{+}$", 0.167: r"$\rm{Cl}^{-}$"}
    Dmax = max(D)
    D0 = [d/Dmax for d in D]
    ax = ax_small
    ax.plot(Z, D0, ".-", label=names[r])
    ax.set_xlabel("z position of molecule [nm]")
    ax.set_ylabel("D/D0")
    ax.set_title("rel. diffusivity (2D model)")
    
# coarser calculation for remaining radii
for r in [0.5, 1., 1.5]:
    N = 2e4
    data = f.get_fields("pugh_diffusivity2D", rMolecule=r, h=4., Nmax=N)
    Z, D = zsorted(data, "D")

    # plot
    ax = ax_big
    ax.plot(Z, D, ".-", label="r=%.3f, N=20k" %r)
    ax.set_xlabel("z position of molecule [nm]")
    ax.set_ylabel("D/D0")
    ax.set_title("rel. diffusivity (2D model)")

ax_big.legend(bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
ax_small.legend(bbox_to_anchor=(1.05, 1.), loc="upper left", borderaxespad=0.,)
ax_small.legend(loc="lower right")

nanopores.savefigs("pugh_diffusivity", folders.FIGDIR)
plt.show()
