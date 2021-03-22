# (c) 2016 Gregor Mitscha-Baude
from nanopores import user_params
import nanopores.models.pughpore as pugh
#from folders import fields
from numpy import linspace
#
#ddata = {2: dict(name="Dpugh", Nmax=1e5, dim=2, r=0.11, h=1.0),
#         3: dict(name="Dpugh", Nmax=2e6, dim=3, r=0.11, h=2.0)}

physp = dict(
    bV = -0.08,
    Qmol = 5.,
    bulkcon = 1000.,
    dnaqsdamp = 0.7353,
)
default = {
    2: dict(physp, dim=2, h=.75, Nmax=1e5, diffusivity="Dpugh2"),
    3: dict(physp, dim=3, h=2., Nmax=6e5, diffusivity="Dpugh2",
                  stokesiter=True)}

dim = 3
params = user_params(default[dim])

X = pugh.tensorgrid(nz=30, nr=4)
X = [[0., 0., t] for t in linspace(-27., 27., 20)]
#result = pugh.F_explicit(X, nproc=2, **params)
result = pugh.F_explicit([None], name="pugh_current_open", **params)
print(result)
print(result["J"])
