''' proof of concept for force interpolation on not-aligned mesh '''

from nanopores import *
from dolfin import *
parameters["allow_extrapolation"] = True

name = "H_geo"
geo_dict = generate_mesh(4e-1, name, dim=2, x0=None)
mesh_file = "%s/%s/mesh/mesh.xml" %(DATADIR,name)

geo1 = geo_from_name(name)
mesh = geo1.mesh
#plot(mesh)
#plot(geo1.subdomains)

print "# Cells:",mesh.num_cells()
print "# Vertices:",mesh.num_vertices()

geo = geo_from_subdomains(mesh, "inside")
plot(geo.subdomains, title="subdomains 1")
mesh = geo.submesh("moleculeinfluid")

geo_i = geo_from_subdomains(mesh, "inside-outside")
plot(geo_i.subdomains, title="subdomains 2")
mesh_i = geo_i.submesh("moleculeinfluid")

DG_i = FunctionSpace(mesh_i, "DG", 0)
f_i = Function(DG_i)
f_i.interpolate(Expression("pow((1-1e8*x[0]),4)"))

V = FunctionSpace(mesh, "CG", 1)
f = Function(V)
f.interpolate(f_i)

plot(mesh_i)
plot(mesh)
plot(sqrt(f_i**2))
plot(f)

interactive()

