# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from nanopores import user_params, user_param
import nanopores.models.pughpore as pugh
from folders import fields

ddata = {2: dict(name="Dpugh", Nmax=1e5, dim=2, r=0.11, h=1.0),
         3: dict(name="Dpugh", Nmax=2e6, dim=3, r=0.11, h=2.0)}

physp = dict(
    bV = -0.08,
    Qmol = 5.,
    bulkcon = 1000.,
    dnaqsdamp = 0.5882,
)
default = {
    2: dict(physp, dim=2, h=.75, Nmax=1e5, diffusivity_data=ddata[2]),
    3: dict(physp, dim=3, h=1.25, Nmax=7e5, diffusivity_data=ddata[3],
                  stokesiter=True)}

dim = user_param(dim=3)
params = user_params(default[dim])
x0 = user_param(x0=[0.,0.,0.])
cache = user_param(cache=True)

result = pugh.F_explicit([x0], name="pugh_vsc_test", cache=cache, **params)

print(result)
print(result["J"])
