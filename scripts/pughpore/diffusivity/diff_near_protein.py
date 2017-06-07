# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import os, sys
from nanopores import fields
from nanopores.models import diffusion_interpolation as diff
from nanopores.models import pughpore as pugh
import dolfin
fields.set_dir_dropbox()

HOME = os.path.expanduser("~")
sys.path.append(HOME + "/Dropbox/nanopores/fenicstools")
from fenicstools import Probes

# set up domain with protein
setup = pugh.Setup(x0=[0.,0.,0.], dim=3, h=8.)
plotter = pugh.Plotter(setup)
meshp = setup.geo.mesh # mesh WITH protein
x0 = np.array(setup.geop.x0)
r0 = setup.geop.rMolecule
rion = 0.11

# load diffusivity on mesh without protein
functions, mesh = fields.get_functions("Dpugh", dim=3, Nmax=2e6)
dist = functions["dist"]
D0 = functions["D"]

# evaluate dist, D on meshp nodes
V = dolfin.FunctionSpace(mesh, "CG", 1)
VV = D0.function_space()

x = meshp.coordinates()
probes = Probes(x.flatten(), V)
probes(dist)
distp = probes.array(0)

probes_ = Probes(x.flatten(), VV)
probes_(D0)
D0p = probes_.array(0)

# initalize D
VVV = dolfin.TensorFunctionSpace(meshp, "CG", 1, shape=(3, 3))
D = dolfin.Function(VVV)
D.interpolate(dolfin.Constant(np.eye(3)))

# first create (N,3,3) array from D0 (N,3)
N = len(D0p)
Da = np.zeros((N, 3, 3))
i3 = np.array([0, 1, 2])
Da[:, i3, i3] = D0p

R = x - x0
r = np.sqrt(np.sum(R**2, 1))
h = r - r0
near = h < distp
#D[~near] =


"""
plotter = Plotter()
mesh2D = plotter.mesh2D
V2D = dolfin.FunctionSpace(mesh2D, "CG", 1)
v2d = dolfin.vertex_to_dof_map(V2D)


probes(dxDx)
y = probes.array(0)
dxDx_2D = dolfin.Function(V2D)
dxDx_2D.vector()[v2d] = y[:]

probes(Dx)
y = probes.array(1)
Dx_2D = dolfin.Function(V2D)
Dx_2D.vector()[v2d] = y[:]




d = D0.compute_vertex_values()
dofs = [VV.sub(i).dofmap().dofs() for i in [0, 1, 2]]
di = [d[dofs[i]] for i in [0, 1, 2]]



dist = dist.compute


def D(x):
    R = np.array(x) - np.array(x0)
    dist_protein = np.sqrt(np.sum(R**2))
    if dist(x) <= dist_protein:
        return D0(x)
    else:
        # treat protein as a plane
        h = dist_protein - r + 1e-100
        Dt = Dt_plane(h, rion)
        Dn = Dn_plane(h, rion)
        Dx = transformation(R, Dn, Dt)
        return Dx
"""