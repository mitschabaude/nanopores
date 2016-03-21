" assess mesh quality of howorka geometry and how it changes with uniform refinement and snapping "

from nanopores import *
from nanopores.geometries.curved import Cylinder, Sphere, Circle
from dolfin import *
from matplotlib import pyplot

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 3.,
z0 = 2.,
ratio = .1,
nref = 1,
)

geo_params = dict(
x0 = [0., 0., nm*z0],
rMolecule = nm*0.5,
lcCenter = nm*0.1, #1,
lcMolecule = nm*0.05, #025,
)

generate_mesh(h, geo_name, optimize=True, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
print geo._bou2phys
#plot(geo.submesh("pore"))
plot_sliced(geo)

# define sphere for molecule
molec = Sphere(R=geo.params["rMolecule"], center=geo.params["x0"])
# define cylinders for inner and outer DNA boundary and side boundary
innerdna = Cylinder(R=geo.params["r0"], L=geo.params["l0"])
outerdna = Cylinder(R=geo.params["r1"], L=geo.params["l0"])
side = Cylinder(R=geo.params["R"], L=2.*geo.params["Rz"])
curved = dict(
    moleculeb = molec.snap,
    innerdnab = innerdna.snap,
    outerdnab = outerdna.snap,
    membranednab = outerdna.snap,
    sideb = side.snap,
    outermembraneb = side.snap,
    )

    
def mesh_quality(mesh, oldmesh=None, ratio=1e-1):
    #vertex = VertexFunction("bool", mesh, False)
    dgncells = CellFunction("size_t", mesh, 0)
    for c in cells(mesh):
        if c.radius_ratio() < ratio:
            dgncells[c] = 1
            #if c.radius_ratio() < 1e-5:
                #print 'Degenerate cell', c.index(), ', radius ratio', c.radius_ratio()
            #for v in vertices(c):
                #vertex[v] = True
                #if c.radius_ratio() < 1e-6:
                #    print '  ', v.point().str()
                
    minrr = MeshQuality.radius_ratio_min_max(mesh)[0]
    print "Minimal radius ratio of mesh:", minrr
    pyplot.figure()
    exec(MeshQuality.radius_ratio_matplotlib_histogram(mesh, 200), locals())
    # plot degenerate cells
    if minrr < ratio:
        submesh = SubMesh(mesh, dgncells, 1)
        title = "degenerate N=%s" %mesh.num_cells()
        #plot(submesh, title=title)
        geo_sub = geo_from_subdomains(submesh,
                    "nanopores.geometries.%s.subdomains" %geo.params["name"], **geo.params)
        plot(geo_sub.boundaries, title="boundaries "+title)
        # find degenerate cells before snapping
        if oldmesh is not None:
            oldmesh = refine(oldmesh)
            oldcells = CellFunction("size_t", oldmesh, 0)
            oldcells.array()[:] = dgncells.array()
            plot(SubMesh(oldmesh, oldcells, 1), "old degenerate cells N=%s" %mesh.num_cells())

    
# mesh quality before refinement
mesh = geo.mesh
print "Number of cells:", mesh.num_cells()
mesh_quality(mesh, ratio=ratio)
#interactive()

for i in range(nref):

    # Mark cells for refinement
    markers = CellFunction("bool", mesh, True)

    # Refine mesh
    mesh = refine(mesh, markers)
    print "Number of cells:", mesh.num_cells()
    geo.adapt(mesh)
    mesh_quality(mesh, ratio=ratio)
    
    # snap curved boundaries
    for boundary, snap in curved.items():
        print "Adapting curved boundary '%s'." % boundary
        geo.snap_to_boundary(boundary, snap)
    mesh_quality(mesh, ratio=ratio)
    
    #areCh = assemble(Constant(1.)*geo.dS("dnab"))
    #print "Area (approx):", areCh
    #print "Error A:", abs(areCh - areC)
 
print "hmin [nm]: ", geo.mesh.hmin()/nm

plot_sliced(geo)
interactive()
