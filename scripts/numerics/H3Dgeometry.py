from nanopores import *
from nanopores.geometries.curved import Cylinder, Sphere, Circle
from dolfin import *

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 1.,
)

geo_params = dict(
x0 = [1., 0., 6.],
rMolecule = 0.5*nm,
lcCenter = 0.1, #5,
lcMolecule = 0.05, #025,
#moleculeblayer = True,
)

# 3D geometry
meshgen_dict = generate_mesh(h, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
plot(geo.submesh("solid"))
plot_sliced(geo)

interactive()
