# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from matplotlib import pyplot as plt
from nanopores import fields, user_params
up = user_params(diamPore=7.)

data = fields.get_fields("diffz_pugh", **up)
z = map(lambda x: x[2], data["x"])
data, z = fields._sorted(data, z)
plt.plot(z, data["D"], "o-")
plt.ylim(0, 1)
