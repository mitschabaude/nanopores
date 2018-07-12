# (c) 2017 Gregor Mitscha-Baude
from nanopores.models.pughpore import Plotter
from nanopores.models.diffusion_interpolation import get_pugh_diffusivity_alt
import dolfin
import os

import nanopores.tools.fields as fields
HOME = os.path.expanduser("~")
DATADIR = os.path.join(HOME, "Dropbox", "nanopores", "fields")
fields.set_dir(DATADIR)
#fields.show()

#params = dict(Nmax=2e5, dim=3, r=2.0779, h=2.0)
# also available (larger mesh):
params = dict(Nmax=5e5, dim=3, r=2.0779, h=1.0)

functions, mesh = get_pugh_diffusivity_alt(**params)
Dfunc = functions["D"]
dis = functions["dist"]
#print "available functions:", functions.keys()
V = dolfin.FunctionSpace(mesh, "CG", 1)

Dx = Dfunc[0]
Dy = Dfunc[1]
Dz = Dfunc[2]
dxDx = dolfin.project(dolfin.grad(Dfunc[0])[0], V)
dyDy = dolfin.project(dolfin.grad(Dfunc[1])[1], V)
dzDz = dolfin.project(dolfin.grad(Dfunc[2])[2], V)

#plotter = Plotter()
#plotter.plot(Dx, title="Dx")
#plotter.plot(dxDx, title="dxDx")
#dolfin.interactive()
