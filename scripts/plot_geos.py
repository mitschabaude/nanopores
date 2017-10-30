# (c) 2017 Gregor Mitscha-Baude
import dolfin
from nanopores.geometries import allpores
#from nanopores import plot_sliced

def plot_sliced(geo, **params):
    tol = 1e-5
    class Back(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return x[1] >= -tol
    back = dolfin.CellFunction("size_t", geo.mesh, 0)
    Back().mark(back, 1)
    submesh = dolfin.SubMesh(geo.mesh, back, 1)
    #plot(submesh)
    bb = geo.mesh.bounding_box_tree()
    subsub = dolfin.CellFunction("size_t", submesh, 0)
    sub = geo.subdomains
    for cell in dolfin.cells(submesh):
        iparent = bb.compute_first_entity_collision(cell.midpoint())
        subsub[cell] = sub[int(iparent)]
    plot_params = dict(title="sliced geometry with subdomains",
                       elevate=-90., **params)
    dolfin.plot(subsub, **plot_params)

geo1 = allpores.get_geo(geoname="wei", subs="solid", h=10., x0=[0.,0.,0.], dim=3)
plot_sliced(geo1, scalarbar=False)
print geo1

geo2 = allpores.get_geo(geoname="alphahem", subs="solid", h=1., x0=[0.,0.,0.], dim=3)
plot_sliced(geo2, scalarbar=False)

dolfin.interactive()