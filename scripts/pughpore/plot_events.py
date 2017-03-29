# -*- coding: utf-8 -*
# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from nanopores import fields
fields.set_dir_dropbox()

from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

drop, t = fields.get("events_pugh_experiment", "drop", "t")
drop, t = np.array(drop), 1e3*np.array(t)
plt.scatter(t, drop, facecolors='#86D466', s=25, marker='s',edgecolors='none',
            label="experiment")

drop, t = fields.get("events3_nobind_new", "a", "t")
drop, t = np.array(drop), 1e3*np.array(t)
plt.scatter(t, drop, facecolors='#027401', s=25, marker='s',edgecolors='none',
            label="simulation")

plt.xscale('log')
plt.xlim(0.02, 2e5)
plt.xlabel(u'τ off [µs]')
plt.ylabel(r'$A/I_0$ (%)')

ax = plt.gca()
ax.set_ylim([0, 40])
ax.set_yticks([0., 10., 20., 30., 40.])
ax.invert_yaxis()
xfmt = FormatStrFormatter("%g")
ax.xaxis.set_major_formatter(xfmt)

plt.legend(loc="upper right", frameon=False)

#plt.tight_layout()
plt.show()