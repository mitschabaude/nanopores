# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import dolfin, os
import nanopores
from nanopores.models.pughpore import Plotter
from nanopores.models.diffusion_interpolation import get_pugh_diffusivity_alt
from nanopores.tools import fields
fields.set_dir(os.path.expanduser("~") + "/Dropbox/nanopores/fields")

params = nanopores.user_params(dim=3, r=2.0779, h=1., Nmax=.5e6)
functions = get_pugh_diffusivity_alt(**params)
D = functions["D"]

print(D([0.,0.,0.]))

#D = dolfin.as_matrix(np.diag([D[i] for i in range(params.dim)]))
plotter = Plotter()
plotter.plot(functions["D"][0], title="Dx")
plotter.plot(functions["D"][1], title="Dz")
plotter.plot(functions["dist"], title="dist")
dolfin.interactive()
