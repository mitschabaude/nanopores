import nanopores.models.pughpore as pughm
from nanopores.tools import fields
import nanopores
import folders

nanopores.add_params(h=1.)
functions, mesh = fields.get_functions("pugh_distance", h=h)
y = functions["y"]
pughm.Plotter().plot(y, title="distance", interactive=True)
