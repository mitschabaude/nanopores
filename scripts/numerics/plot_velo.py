# (c) 2016 Gregor Mitscha-Baude
import matplotlib.pyplot as plt
import numpy as np
import nanopores.tools.fields as fields
import os
import nanopores

DATADIR = os.path.join(nanopores.HOME, "Dropbox", "nanopores", "fields")
FIGDIR = os.path.join(nanopores.HOME, "Dropbox", "nanopores", "figures")
FIGDIR = os.path.join(nanopores.HOME, "papers", "pnps-numerics", "figures")
fields.set_dir(DATADIR)
fields.update()
#data = fields.get_fields("howorka_velo1", Nmax=1e4)
data = fields.get_fields("howorka_velo3D_2")

#z = [x[2] for x in data["x"]]
z = [x[0] for x in data["x"]]
data, z = fields._sorted(data, z)

## modify
i = 0
rad = lambda x: sum(t**2 for t in x)
#rad = lambda x: x[2]
#v0 = [rad(x) for x in data["v0"]]
#v1 = [rad(x) for x in data["v1"]]
z.pop(-1)
data["v0"].pop(-1)
data["v1"].pop(-1)

plt.plot(z, [x[0] for x in data["v0"]], "o-", color="#000099", label=r"$v_r$ (exact)")
plt.plot(z, [x[0] for x in data["v1"]], ".--", color="#0099ff", label=r"$v_r$ (approx.)")
plt.xlabel("radial position of molecule [nm]") ## modify
plt.ylabel(r"velocity [m/s]")
plt.legend(loc="best") #(bbox_to_anchor=(0.5, 1.), loc="upper left")
#plt.plot(z, np.abs(np.array(v0)-np.array(v1)))
fig = plt.gcf()
fig.set_size_inches((4,3))

plt.figure()

plt.plot(z, [x[2] for x in data["v0"]], "v-", color="#009900", label=r"$v_z$ (exact)")
plt.plot(z, [x[2] for x in data["v1"]], ".--", color="#99ff00", label=r"$v_z$ (approx.)")
plt.xlabel("radial position of molecule [nm]") ## modify
plt.ylabel("velocity [m/s]")
plt.legend(loc="best") #bbox_to_anchor=(0., 1.), loc="upper left")
#plt.plot(z, np.abs(np.array(v0)-np.array(v1)))
fig = plt.gcf()
fig.set_size_inches((4,3))

nanopores.savefigs("velo", FIGDIR)
#plt.show()