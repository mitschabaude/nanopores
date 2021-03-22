# (c) 2016 Gregor Mitscha-Baude
import diffusion
import nanopores.tools.solvers as solvers
import nanopores.models.pughpore as pugh
import numpy as np
from folders import fields
import nanopores as nano
import dolfin

def tolist(array):
    return [list(a) for a in array]

up = nano.user_params(h=1., Nmax=2e5, H=10., R=3.1)
params = dict(
    H = up.H,
    R = up.R,
    center_z_at_x0 = True, #False, #True,
    dim = 3,
    h = up.h,
    Nmax = up.Nmax,
    rMolecule = 2.0779,
    bulkbc = True,
)

@solvers.cache_forcefield("pugh_diff3D", {})
def D_tensor(X, **params):
    x0 = X[0]
    print()
    print()
    print("MOLECULE: ",x0)
    setup = pugh.Setup(x0=x0, **params)
    nano.plot_sliced(setup.geo)
    #dolfin.plot(setup.geo.subdomains, elevate=-45., key="subdomains")
    #dolfin.plot(setup.geo.subdomains, key="subdomains")
    setup.physp["bulkbc"] = params["bulkbc"]
    D = diffusion.diffusivity_tensor(setup)
    return dict(D=[tolist(D)])

def calculate_diff_field(params):
    X = fields.get_entry("pughx", "x", k=3)
    return D_tensor(X, nproc=1, **params)

def calculate_1D_profile(params, N=100, H=50.):
    Z = np.linspace(-H, H, N)
    X = [[0.,0.,z] for z in Z]
    params2D = dict(params, h=1., Nmax=1e5, H=150., R=75.)
    diffusion.calculate_diffusivity2D(X, nproc=5, name="pugh_diff2D", **params2D)

def calculate_diff_plot(params):
    # calculate 2D reference value
    params2D = dict(params, h=.1, Nmax=1e5, H=16.)
    diffusion.calculate_diffusivity2D([[0.,0.,0.]],
         name="pugh_diff2D_test", **params2D)

    # create points for 3D
    l0 = pugh.pughpore.params["l3"]
    r = params["rMolecule"]
    eps = 1e-2
    R = l0/2. - r - eps
    X = [[t, 0., 0.] for t in np.linspace(0, R, 20)]

    # calculate
    D_tensor(X, name="pugh_diff3D_test", nproc=2, **params)

# TODO: next time, use geometry with x=0 corresponding to the wall!!
def calculate_D_outside(params):
    # create points for 3D
    eps = 1e-2
    l0 = pugh.pughpore.params["l0"]
    r = params["rMolecule"]
    eps = 1e-2
    R0 = l0/2. + r + eps
    R1 = (l0/2. + params["R"])/2.
    X = [[t, 0., 0.] for t in np.linspace(R0, R1, 40)]
    D_tensor(X, name="pugh_diff3D_test", nproc=5, **params)

def calculate_D_inside(params):
    # create points for 3D
    eps = 1e-2
    l3 = pugh.pughpore.params["l3"]
    r = params["rMolecule"]
    eps = 1e-2
    R0 = 0.
    R1 = l3/2. - r - eps
    X = [[t, 0., -1.] for t in np.linspace(R0, R1, 30)]
    D_tensor(X, name="pugh_diff3D_cross", nproc=5,
             dnaqsdamp=0.0, lcMolecule=.4, cheapest=True, **params)

if __name__ == "__main__":
    pass
    #D_tensor([[params["R"]/2. + 9./2.,0.,0.]], cache=False,
    #         name="pugh_diff3D_test", nproc=1, **params)
    #calculate_diff_plot(params)
#    X = [[0.,0.,40.]]
#    params2D = dict(params, h=1., Nmax=1e5, H=150., R=75.)
#    diffusion.calculate_diffusivity2D(X, nproc=1, cache=False,
#                                      name="pugh_diff2D", **params2D)

    #calculate_1D_profile(params)
    #calculate_D_outside(params)
#    D_tensor([[0.5, 0., -1.]], cache=False, dnaqsdamp=0.0, lcMolecule=.4,
#             cheapest=True,
#             name="pugh_diff3D_cross", nproc=1, **params)
#    dolfin.interactive()
    calculate_D_inside(params)