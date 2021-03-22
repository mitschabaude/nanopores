# (c) 2017 Gregor Mitscha-Baude
from itertools import product
import numpy as np
import dolfin
from nanopores.tools import fields

def Dt_plane(h, r):
    x = r/h
    return 1. - 9./16.*x + 1./8.*x**3 - 45./256.*x**4 - 1/16.*x**5

sinh = np.sinh
acosh = np.arccosh
def Dn_plane(l, r, N=100):
    alpha = acosh(l/r)
    s = 0.
    for n in range(1, N):
        n = float(n)
        K = n*(n+1)/(2*n-1)/(2*n+3)
        s += K*((2*sinh((2*n+1)*alpha)+(2*n+1)*sinh(2*alpha))/(4*(sinh((n+.5)*alpha))**2-(2*n+1)**2*(sinh(alpha))**2) - 1)
    return 1./((4./3.)*sinh(alpha)*s)

def set_D_with_protein(setup):
    meshp = setup.geo.mesh # mesh WITH protein
    x0 = np.array(setup.geop.x0)
    dim = setup.phys.dim
    x0 = x0 if dim==3 else x0[::2]
    r0 = setup.geop.rMolecule
    rion = 0.11

    # load diffusivity on mesh without protein
    functions, mesh = fields.get_functions(**setup.solverp.diffusivity_data)
    #dist = functions["dist"]
    D0 = functions["D"]

    # evaluate dist, D on meshp nodes
    #Vp = dolfin.FunctionSpace(meshp, "CG", 1)
    VVp = dolfin.VectorFunctionSpace(meshp, "CG", 1)

    #distp_ = dolfin.interpolate(dist, Vp)
    D0p_ = dolfin.interpolate(D0, VVp)

    #distp = distp_.vector()[dolfin.vertex_to_dof_map(Vp)]
    D0p = D0p_.vector()[dolfin.vertex_to_dof_map(VVp)]
    D0p = np.column_stack([D0p[i::dim] for i in range(dim)])

    x = meshp.coordinates()

    # first create (N,3,3) array from D0 (N,3)
    N = len(D0p)
    Da = np.zeros((N, dim, dim))
    i3 = np.array(list(range(dim)))

    Da[:, i3, i3] = D0p

    # modify (N,3,3) array to include protein interaction
    R = x - x0
    r = np.sqrt(np.sum(R**2, 1))
    overlap = r < rion + r0
    #near = ~overlap & (r - r0 < distp)
    h = np.maximum(r - r0, rion)
    eps = 1e-2
    D00 = setup.phys.D

    Dt = np.zeros_like(r)
    Dn = np.zeros_like(r)

    Dt[overlap] = eps
    Dn[overlap] = eps
    Dt[~overlap] = Dt_plane(h[~overlap], rion)
    Dn[~overlap] = Dn_plane(h[~overlap], rion, N=20)

    # D(R) = Dn(h) RR^T + Dt(h) (I - RR^T) where R is normalized
    R0 = R/(r[:, None] + 1e-100)
    RR = (R0[:, :, None] * R0[:, None, :])
    I = np.zeros((N, dim, dim))
    I[:, i3, i3] = 1.
    Dpp = D00*(Dn[:, None, None] * RR + Dt[:, None, None] * (I - RR))
    near = Dpp[:, dim-1, dim-1] < Da[:, dim-1, dim-1]
    Da[near] = Dpp[near]

    # assign final result to dolfin P1 TensorFunction
    VVV = dolfin.TensorFunctionSpace(meshp, "CG", 1, shape=(dim, dim))
    D = dolfin.Function(VVV)
    v2d = dolfin.vertex_to_dof_map(VVV)
    Dv = D.vector()
    for k, (i,j) in enumerate(product(i3, i3)):
        Dv[np.ascontiguousarray(v2d[k::dim**2])] = np.ascontiguousarray(Da[:, i, j])

    setup.phys.update(Dp=D, Dm=D)
    return D