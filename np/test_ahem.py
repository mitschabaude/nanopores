from nanopores import *
from dolfin import *

geo_params = dict(
    l4 = 10.,
    R = 20.
)

t = Timer("meshing")
meshdict = generate_mesh(6., "aHem", **geo_params)

print "Mesh generation time:",t.stop()
print "Mesh file:",meshdict["fid_xml"]
print "Mesh metadata:"
for item in meshdict["meta"].items():
    print "'%s' : %s" %item
print 

t = Timer("reading geometry")
geo = geo_from_xml("aHem")

print "Geo generation time:",t.stop()
print "Geo params:", geo.params
print "Geo physical domains:", geo._physical_domain
print "Geo physical boundaries:", geo._physical_boundary

plot(geo.submesh("ahem"))
plot(geo.submesh("fluid_bulk_top"))
plot(geo.submesh("fluid_center"))
interactive()
