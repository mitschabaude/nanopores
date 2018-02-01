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
    #dolfin.plot(subsub, **plot_params)
    return subsub

params = dict(geoname="wei", subs="solid", h=5., x0=[0.,0.,0.], dim=3, rMolecule=1.25)
zreceptor = 0.95
zprotein = 0.90
rrec = 0.5
distrec = 4. - params["rMolecule"] - rrec
pore = allpores.get_pore(**params)
ztop = pore.protein.zmax()[1]
zbot = pore.protein.zmin()[1]
zrec = zbot + rrec + (ztop - zbot - 2.*rrec)*zreceptor
xrec = pore.radius_at(zrec) - distrec
zprot = zbot + params["rMolecule"] + (ztop - zbot - 2.*params["rMolecule"])*zprotein
params["receptor"] = [xrec, 0., zrec]
params["rReceptor"] = rrec
params["x0"] = [8., 3., zprot]
geo1 = allpores.get_geo(**params)
domains1 = plot_sliced(geo1, scalarbar=False)
print geo1

geo2 = allpores.get_geo(geoname="alphahem", subs="solid", h=.5, x0=None, dim=3)
domains2 = plot_sliced(geo2, scalarbar=False)

from nanopores.dirnames import DROPBOX
file1 = dolfin.File(DROPBOX + "/geo_wei.pvd")
file2 = dolfin.File(DROPBOX + "/geo_ahem.pvd")

#file1 << geo1.mesh
file1 << geo1.subdomains

#file2 << geo2.mesh
file2 << geo2.subdomains
#dolfin.interactive()