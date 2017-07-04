# (c) 2017 Gregor Mitscha-Baude
from nanopores import geo_from_name, generate_mesh, user_params, import_vars

geo_name = "W_2D_geo"
geo_dict = import_vars("nanopores.geometries.%s.params_geo" %geo_name)
physical_dict = import_vars("nanopores.physics.params_physical")
default_dict = dict(geo_dict = geo_dict, physical_dict = physical_dict)

nm = geo_dict["nm"]
lsam = geo_dict["lsam"]
params = dict(
    x0 = [0, 0, 0],
    outerfrac = 0.2,
    lsam = lsam,
    r0 = 13*nm-lsam,
    angle = 40,
    sam = None,
)

generate_mesh(1., geo_name, **params)
geo = geo_from_name(geo_name, **params)

geo.plot_boundaries()
geo.plot_subdomains(interactive=True)

