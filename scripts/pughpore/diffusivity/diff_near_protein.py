# (c) 2017 Gregor Mitscha-Baude
import numpy as np
#import os, sys
from itertools import product
from nanopores import fields
from nanopores.models import diffusion_interpolation as diff
from nanopores.models import pughpore as pugh
import dolfin
fields.set_dir_dropbox()

#HOME = os.path.expanduser("~")
#sys.path.append(HOME + "/Dropbox/nanopores/fenicstools")
#from fenicstools import Probes

def set_D_with_protein(setup):
    meshp = setup.geo.mesh # mesh WITH protein
    x0 = np.array(setup.geop.x0)
    dim = setup.phys.dim
    x0 = x0 if dim==3 else x0[::2]
    r0 = setup.geop.rMolecule
    rion = 0.11

    # load diffusivity on mesh without protein
    functions, mesh = fields.get_functions(**setup.solverp.diffusivity_data)
    dist = functions["dist"]
    D0 = functions["D"]

    # evaluate dist, D on meshp nodes
    Vp = dolfin.FunctionSpace(meshp, "CG", 1)
    VVp = dolfin.VectorFunctionSpace(meshp, "CG", 1)

    distp_ = dolfin.interpolate(dist, Vp)
    D0p_ = dolfin.interpolate(D0, VVp)

    distp = distp_.vector()[dolfin.vertex_to_dof_map(Vp)]
    D0p = D0p_.vector()[dolfin.vertex_to_dof_map(VVp)]
    D0p = np.column_stack([D0p[i::dim] for i in range(dim)])

    x = meshp.coordinates()
    #    probes = Probes(x.flatten(), V)
    #    probes(dist)
    #    distp = probes.array(0)
    #
    #    probes_ = Probes(x.flatten(), VV)
    #    probes_(D0)
    #    D0p = probes_.array(0)

    # first create (N,3,3) array from D0 (N,3)
    N = len(D0p)
    Da = np.zeros((N, dim, dim))
    i3 = np.array(list(range(dim)))

    Da[:, i3, i3] = D0p

    # modify (N,3,3) array to include protein interaction
    R = x - x0
    r = np.sqrt(np.sum(R**2, 1))
    overlap = r < rion + r0
    near = ~overlap & (r - r0 < distp)
    h = np.maximum(r[near] - r0, rion)
    eps = 1e-2
    D00 = setup.phys.D

    Dt = np.zeros_like(r)
    Dn = np.zeros_like(r)

    Dt[overlap] = eps
    Dn[overlap] = eps
    Dt[near] = diff.Dt_plane(h, rion)
    Dn[near] = diff.Dn_plane(h, rion, N=20)

    # D(R) = Dn(h) RR^T + Dt(h) (I - RR^T) where R is normalized
    R0 = R/(r[:, None] + 1e-100)
    RR = (R0[:, :, None] * R0[:, None, :])
    I = np.zeros((N, dim, dim))
    I[:, i3, i3] = 1.
    Dpp = Dn[:, None, None] * RR + Dt[:, None, None] * (I - RR)
    Da[overlap | near] = D00*Dpp[overlap | near]

    # assign final result to dolfin P1 TensorFunction
    VVV = dolfin.TensorFunctionSpace(meshp, "CG", 1, shape=(dim, dim))
    D = dolfin.Function(VVV)
    v2d = dolfin.vertex_to_dof_map(VVV)
    Dv = D.vector()
    for k, (i,j) in enumerate(product(i3, i3)):
        Dv[np.ascontiguousarray(v2d[k::dim**2])] = np.ascontiguousarray(Da[:, i, j])

    setup.phys.update(Dp=D, Dm=D)
    #embed()
    return D


# set up domain with protein
ddata = dict(name="Dpugh", r=0.11, dim=2)
setup = pugh.Setup(x0=[0., 0., 10.], dim=2, h=1., Nmax=1e5,
                   diamPore=6.*np.sqrt(np.pi)/2.,
                   cheapest=True, diffusivity_data=ddata)
plotter = pugh.Plotter(setup)
pugh.prerefine(setup, True)

print("getting D")
D = set_D_with_protein(setup)
print("plotting")
plotter.plot(D[0, 0], title="Dxx")
dolfin.interactive()