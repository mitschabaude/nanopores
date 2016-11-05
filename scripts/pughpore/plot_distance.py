from dolfin import *
import nanopores.models.pughpore as pughm
from nanopores.tools import fields

fields.set_dir("/home/benjamin/Desktop/data/")
functions, mesh = fields.get_functions("pugh_distance")
dis = functions["pugh_distance"]
setup=pughm.Setup()
plotter=pughm.Plotter(setup)
plotter.plot(dis,title="Distance")
interactive()
