"""
python tools for creating axisymmetric geometries automatically from their half-plane.
"""
from .. import py4gmsh
from . import box
from .box import BoxCollection, entities_to_gmsh, to_mesh
import dolfin

__all__ = ["rotate_z"]

class CylinderCollection(BoxCollection):
    """complement BoxCollection for creation of 3D cylinders.
    intended usage:
        
    >> box = Box([0., 0.], [1., 1.])
    >> subbox = Box([1./3, 1./3], [2./3, 2./3])
    >> domain = box - subbox
    >> ... some further operations specifying the domain, e.g. boundaries ...
    >> domain3D = rotate_z(domain) # rotate around second axis
    >> geo = domain3D.create_geometry(lc=.1)
    
    rotate_z() takes a BoxCollection and returns a CylinderCollection.
    CylinderCollection will simply inherit everything from the given
    BoxCollection, but will act slightly differently when create_geometry
    is called; it will create a 3D gmsh geometry by invoking the
    py4gmsh.Extrude() function at some point.
    """
        
    def create_geometry(self, lc=.5):
        self.compute_entities()
        gmsh_entities = entities_to_gmsh(self.entities, self.indexsets, lc=lc)
        rotate_surfs(gmsh_entities, self.boundaries, 1)
        rotate_surfs(gmsh_entities, self.subdomains, 2)
        self.geo = to_mesh()
        self.geo.params = self.params
        if hasattr(self, "synonymes"):
            self.geo.import_synonymes(self.synonymes)
        return self.geo
        
Surf = {0 : "Point", 1 : "Line", 2 : "Surface"}
PhysSurf = {
    1 : lambda vol, name: py4gmsh.PhysicalSurface(vol, name, 3),
    2 : lambda vol, name: py4gmsh.PhysicalVolume(vol, name, 3)
} 
        
def rotate_surfs(gmsh_entities, subdomains=[], k=2):
    # this rotates surfaces or lines with Extrude around z axis
    rot_axis = [0.0, 0.0, 1.0]
    point_on_axis = [0.0, 0.0, 0.0]
    angle = 'Pi/2'
    n_rotations = 4
    surfs = gmsh_entities[k]
    vols = [[] for s in surfs]
    
    for i, surf in enumerate(surfs):
        for j in range(n_rotations):
            name = py4gmsh.Extrude("%s{%s}" % (Surf[k], surf), rotation_axis=rot_axis,
                point_on_axis=point_on_axis, angle=angle)
            vols[i].append(name + '[1]')
            surf = name + '[0]'
            
    for sub in subdomains:
        physvol = [vol for i in sub.indexset for vol in vols[i]]
        PhysSurf[k](physvol, sub.name)
    return vols

def rotate_z(domain):
    """take BoxCollection and return equivalent CylinderCollection by
    rotating about the second axis. thus, transform coordinates of
    points like (x, z) --> (x, 0, z)."""
    return rotate(domain, d=1)

def rotate(domain, d=2):
    # TODO: this shouldn't mess with original domain
    # create empty CylinderCollection
    cyl = CylinderCollection()
    
    # copy everything from domain #FIXME
    cyl.__dict__ = dict(domain.__dict__)

    # add csg to make it compatible with Box
    cyl.csg = domain._csg()
    
    # add dummy y dimension to make it 3D
    for box in cyl.boxes + cyl.facets:
        add_dim(box, d)
    return cyl
        
def add_dim(box, d):
    # add dummy dimension to box
    box.intervals.insert(d, (0., 0.))
    box.a = tuple(x for x,y in box.intervals)
    box.b = tuple(y for x,y in box.intervals)
    box.dim = box.dim + 1

