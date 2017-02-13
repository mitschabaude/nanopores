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
    dnaqsdamp = .5
)
default = {
    2: dict(physp, dim=2, h=.75, Nmax=1e5, diffusivity_data=ddata[2]),
    3: dict(physp, dim=3, h=1.5, Nmax=7e5, diffusivity_data=ddata[3],
                  stokesiter=True)}
                  
dim = 2
params = user_params(default[dim])

z = [-26., 1., 26.]
X = [[0.,0.,t] for t in z]

result = pugh.F_explicit(X, nproc=4, name="pugh_test", cache=False, **params)
print result
print result["J"]


