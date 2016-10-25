# (c) 2016 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
from diffusion import calculate_diffusivity2D

def zsorted(data, field):
    z = [x[2] for x in data["x"]]
    J = data[field]
    I = sorted(range(len(z)), key=lambda k: z[k])
    z1 = [z[i] for i in I]
    J1 = [J[i] for i in I]
    return z1, J1
    
# points
H = 42.
Z = np.linspace(-H, H, 48)
X = [[0.,0.,z] for z in Z]

# get data
data = calculate_diffusivity2D(X, nproc=6)
Z, D = zsorted(data, "D")

# plot
fig, ax = plt.subplots(figsize=(5, 4), num="diffusivity")
ax.plot(Z, D, "s-", label="2D")
ax.set_xlabel("z position of molecule [nm]")
ax.set_ylabel("D/D0")
ax.set_title("rel. diffusivity (2D model)")
plt.show()
