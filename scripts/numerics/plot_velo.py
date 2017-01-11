# (c) 2016 Gregor Mitscha-Baude
import matplotlib.pyplot as plt
import numpy as np
import nanopores.tools.fields as fields
import os
import nanopores

DATADIR = os.path.join(nanopores.HOME, "Dropbox", "nanopores", "fields")
FIGDIR = os.path.join(nanopores.HOME, "Dropbox", "nanopores", "figures")
fields.set_dir(DATADIR)
data = fields.get_fields("howorka_velo1", Nmax=1e4)

z = [x[2] for x in data["x"]]
data, z = fields._sorted(data, z)

v0 = [x[1] for x in data["v0"]]
v1 = [x[1] for x in data["v1"]]

plt.plot(z, v1, "o-", label="v (exact)")
plt.plot(z, v0, ".--r", label="v (linear approximation)")
plt.xlabel("z position of molecule [nm]")
plt.ylabel("molecule velocity [m/s]")
plt.legend(loc="upper left")
#plt.plot(z, np.abs(np.array(v0)-np.array(v1)))
fig = plt.gcf()
fig.set_size_inches((5,4))
nanopores.savefigs("howorka_velo", FIGDIR)
#plt.show()