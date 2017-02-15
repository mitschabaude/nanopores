# (c) 2016 Gregor Mitscha-Baude
from nanopores import user_params
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

dim = 3
params = user_params(default[dim])

X = pugh.tensorgrid(nz=30, nr=4)
result = pugh.F_explicit(X, nproc=3, name="pugh_vsc", **params)

print result
print result["J"]

