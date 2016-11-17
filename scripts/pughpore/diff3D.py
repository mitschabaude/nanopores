# (c) 2016 Gregor Mitscha-Baude
import diffusion
import nanopores.tools.solvers as solvers
import nanopores.models.pughpore as pugh
import numpy as np
from matplotlib import pyplot as plt
from itertools import product
#import folders
import nanopores as nano
#default = diffusion.default

def tolist(array):
    return [list(a) for a in array]

up = nano.user_params(h=.4, Nmax=2e5, H=8.)
params = dict(
    H = up.H,
    R = 4.,
    center_z_at_x0 = True,
    dim = 3,
    h = up.h,
    Nmax = up.Nmax,
    rMolecule = 0.11,
    bulkbc = True,
)

@solvers.cache_forcefield("pugh_diff3D_test", {})
def D_tensor(X, **params):
    x0 = X[0]
    setup = pugh.Setup(x0=x0, **params)
    setup.physp["bulkbc"] = params["bulkbc"]
    D = diffusion.diffusivity_tensor(setup)
    return dict(D=[tolist(D)])

# calculate 2D reference value
params2D = dict(params, h=.1, Nmax=1e5, H=16.)
D2D = diffusion.calculate_diffusivity2D([[0.,0.,0.]],
     name="pugh_diff2D_test", **params2D)

# create points for 3D
l0 = pugh.pughpore.params["l3"]
r = params["rMolecule"]
eps = 1e-2
R = l0/2. - r - eps
X = [[t, 0., 0.] for t in np.linspace(0, R, 20)]

# calculate
data = D_tensor(X, nproc=3, **params)

