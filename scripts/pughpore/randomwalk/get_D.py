import nanopores.geometries.pughpore as pughpore
import nanopores as nano
from nanopores.tools import fields
fields.set_dir("/home/benjamin/Dropbox/Paper Howorka/data/fields")

def zsorted(data, field):
    z = [x[2] for x in data["x"]]
    J = data[field]
    I = sorted(range(len(z)), key=lambda k: z[k])
    z1 = [z[i] for i in I]
    J1 = [J[i] for i in I]
    return z1, J1

geop = nano.Params(pughpore.params)
rMolecule = geop.rMolecule
N = 2e4

data = fields.get_fields("pugh_diffusivity2D", rMolecule=rMolecule, h=4., Nmax=N)
Z, D = zsorted(data, "D")
