'''
python class to handle n-dimensional rectangular boxes and unions of them

note: very naive implementation of union, suitable for a couple 100 boxes
'''

# TODO: collaps pathological boxes,
#       reduce number of boxes in disjoint union when there is no intersection

from itertools import izip, product, combinations
import numpy as np
from .. import py4gmsh
import dolfin

class Box(object):
    
    def __init__(self, a=None, b=None, intervals=None):
        if intervals is not None:
            self.intervals = intervals
            self.dim = len(intervals)
        else:
            # a, b should be tuples
            assert len(a) == len(b)
            a, b = _sort(a, b)
            a = tuple(Float(x) for x in a)
            b = tuple(Float(x) for x in b)
            self.intervals = zip(a, b)
            self.dim = len(a)

    def __str__(self):
        return "Box(%s)" %("x".join(["[%s, %s]" %i for i in self.intervals]),)
        
    def __cmp__(self, other):
        return cmp(self.intervals, other.intervals)
        
    def __iter__(self):
        return iter(self.intervals)
        
    def __getitem__(self, i):
        return self.intervals[i]
        
    #def __or__(self, other):
    #    # return union of two Boxes as a BoxCollection
    #    assert isinstance(other, Box)
    #    return multi_box_union(self, other)
    
    @staticmethod    
    def union(*boxes):
        return multi_box_union(*boxes)

        
      
class BoxCollection(object):
    " collection of disjoint boxes, and their vertices, edges, facets "
    
    def __init__(self, *boxes):
        self.boxes = list(boxes)
        #if len(self.boxes):
        #    self.initialize()

    def __str__(self):
        return "Union{\n  %s\n}\n" %(",\n  ".join(str(box) for box in self),)
        
    def __iter__(self):
        return iter(self.boxes)
        
    def __getitem__(self, i):
        return self.boxes[i]    
        
    def append(self, item):
        self.boxes.append(item)
        
    def extend(self, seq):
        self.boxes.extend(seq)
        
    def initialize(self):
        self.boxes = multi_box_union(*(self.boxes)).boxes
        
    def plot(self, clscale=.5, sub=False):
        entities_to_gmsh(self.entities)
        mesh = to_mesh(clscale)
        dolfin.plot(mesh, interactive=True)
        if sub:
            pass
        
    def addboxes(self, *newboxes):
        new = multi_box_union(*(self.boxes + list(newboxes)))
        self.nodes = new.nodes
        self.entities = new.entities
        self.esets = new.esets
        self.boxes = new.boxes
        
    #def __or__(self, other):
    #    if isinstance(other, Box):
    #        return multi_box_union(other, *(self.boxes))
    #        
    #    if isinstance(other, BoxCollection):
    #        return multi_box_union(*(self.boxes + other.boxes))
        
        
class deprecatedPoint(object):
    " n-dimensional point "
    
    def __init__(self, a):
        self.tup = tuple(Float(x) for x in a)
    
    def __cmp__(self, other):
        " test for equality up to tolerance "
        return cmp(self.tup, other.tup)
        
    def __eq__(self, other):
        return cmp(self, other) == 0

    def __ne__(self, other):
        return cmp(self, other) != 0
        
    def __len__(self):
        return len(self.tup)
        
    def __getitem__(self, i):
        return self.tup[i]
        
    def __iter__(self):
        return iter(self.tup)
       

class Float(float):
    tol = 1e-10
    
    def __cmp__(self, other):
        return _cmp(self, other, self.tol)
        
    def __eq__(self, other):
        if isinstance(other, Float):
            return cmp(self, other) == 0
        else:
            return False

    def __ne__(self, other):
        return not self == other
    
        
def _sort(a, b):
    a_ = tuple(map(min, izip(a,b)))
    b_ = tuple(map(max, izip(a,b)))
    return a_, b_
    
    
def _cmp(s, t, tol):
    " test floats for approximate equality "
    if s > t+tol:
        return 1
    elif s < t-tol:
        return -1
    else:
        return 0
               
def _sorted(I, J):
    return (min(I,J), max(I,J))       
        
def _unique(seq):
    # version of list(set(seq)) but using list membership test
    # (weaker than set membership test, for approximate Float equality)
    unique = []
    [unique.append(x) for x in seq if not unique.count(x)]
    return unique
        
def multi_interval_union(*intvs):
    " return disjoint subintervals with same union as intervals, plus parent information "
    # input: list of intervals
    # output: for lists:
    #   1) sorted list of N (unique) endpoints
    #   2) equivalent to 4) for the endpoints
    #   3) corresponding list of N-1 subintervals
    #   4) list of N-1 sets denoting to which of the original intervals a subinterval belongs
    # (number of intervals is expected to be smallish, like < 10000)
    # interval := tuple (a,b) with a <= b, a,b are Floats
    # in the case of a == b the interval will be ignored in the final output
    points = sorted(_unique([x for intv in intvs for x in intv]))
    subs = zip(points[:-1], points[1:])
    psets = [set() for i in points]
    ssets = [set() for i in subs]
    for i, intv in enumerate(intvs):
        a, b = intv
        for j in range(points.index(a), points.index(b)):
            psets[j].add(i)
            ssets[j].add(i)
        psets[points.index(b)].add(i)
    #print points, psets, subs, ssets    
    return points, psets, subs, ssets
    
        
def _intersection(A, B):
    ints = [(max(a,c), min(b,d)) for (a,b),(c,d) in izip(A, B)]
    return Box(intervals=ints)

def _disjoint(A, B):
    " test whether two boxes are disjoint "
    return any((a > d) | (c > b) for (a,b),(c,d) in izip(A, B))
    
    
def multi_box_union(*boxes):
    " return union of boxes as collection of disjoint boxes "
    
    # get dimension; all boxes must have the same
    dim = boxes[0].dim
    assert all([dim == box.dim for box in boxes])
    
    # get list of disjoint intervals for every dimension
    intvs = []
    isets = []
    nodes = []
    nsets = []
    for ilist in izip(*boxes):
        node, nset, intv, iset = multi_interval_union(*ilist)
        nodes.append(node)
        nsets.append(nset)
        intvs.append(intv)
        isets.append(iset)
        
    D = range(dim) # [0,...,dim-1]
    D1 = range(dim+1) # [0,..,dim]
    entities = [[] for k in D1] # entity := d-tuple of tuples/floats <=> box
    esets = [[] for k in D1] # eset := {indices of the parent boxes of an entity}
    
    # collect all admissible k-entities for k = 1,...,d
    # 0-entity = vertex, 1-entity = edge, ...
    for k in D1:
        for tup in combinations(D, k):
            k_entities = product(*[(intvs[i] if i in tup else nodes[i]) for i in D])
            k_esets1D = product(*[(isets[i] if i in tup else nsets[i]) for i in D])
            for entity, eset1D in izip(k_entities, k_esets1D):
                eset = set.intersection(*eset1D)
                if eset:
                    entities[k].append(entity)
                    esets[k].append(eset)
        #print "k=%s:"%k, entities[k]
        #print
                    
    collection = BoxCollection()
    collection.nodes = nodes # TODO: needed?
    collection.entities = entities
    collection.esets = esets
    collection.boxes = [Box(intervals=A) for A in entities[dim]]
            
    return collection
        
def _facets(box): # get facets from d-dimensional box, d >= 1
    facets = []
    d = len(box)
    for i, x in enumerate(box): # iterate through intervals/floats
        if isinstance(x, tuple):
            facets.append(tuple((x[0] if i==j else box[j]) for j in range(d)))
            facets.append(tuple((x[1] if i==j else box[j]) for j in range(d)))
    return tuple(facets)
    
def _orientations(k): # facet orientations
    return ["-","","","-","-",""][:(2*(k+1))] # valid for 0-, 1- and 2-dim. facets
    
FacetLoop = {
    0: lambda x:x,
    1: py4gmsh.LineLoop,
    2: py4gmsh.SurfaceLoop
    }
    
Entity = {
    1: lambda x: py4gmsh.Line(*x),
    2: py4gmsh.PlaneSurface,
    3: py4gmsh.Volume
    }

def entities_to_gmsh(entities, lc=.5):

    global gmsh_entities
    global dim
    
    dim = len(entities)-1
    gmsh_entities = [[None for e in k] for k in entities]
    
    # a shortcut
    _gmsh = lambda entity, k: gmsh_entities[k][entities[k].index(entity)]
    
    # add points
    for i, e in enumerate(entities[0]):
        # py4gmsh.Point expects three values x,y,z
        e = e + tuple(0. for i in range(3 - dim))
        # create point
        gmsh_entities[0][i] = py4gmsh.Point(e, lc)
        
    # add entities of dimension > 0
    for k in range(1, dim+1):
        for i, e in enumerate(entities[k]):
            # get surrounding facets (x,y,z in this order)
            facets = _facets(e)
            
            # find corresponding gmsh facets
            facets = [_gmsh(l, k-1) for l in facets]
            
            # create facet loop with correct orientations
            orient = _orientations(k-1)
            loop = FacetLoop[k-1]([o+s for o, s in izip(orient, facets)])
            
            # add entity
            gmsh_entities[k][i] = Entity[k](loop)
            
            
def subdomains_to_gmsh(esets):
    # call only after entities_to_gmsh
    sets = esets[dim]
    vols = list(set.intersection(*sets))
    #for i in vols:
    #    for vol in gmsh_entities[
    
    
def to_mesh(clscale, pid=""):
    code = py4gmsh.get_code()

    import os, subprocess
    inputfile = "input%s.geo" %pid
    outfile = "out%s.msh" %pid
    meshfile = "mesh%s.xml" %pid

    # create path/to/nanoporesdata/gid/mesh if not already there
    meshdir = "/tmp/nanopores"
    if not os.path.exists(meshdir):
        os.makedirs(meshdir)

    fid_dict = {"fid_geo": os.path.join(meshdir, inputfile),
                "fid_msh": os.path.join(meshdir, outfile)}

    # save code to .geo file
    fobj = open(fid_dict["fid_geo"], "w")
    fobj.write(code)
    fobj.close()

    # after writing the geo file, call gmsh
    gmsh_out = subprocess.call(["gmsh", "-3", "-v", "1","-clscale", "%f" %clscale,
                     fid_dict["fid_geo"], "-o", fid_dict["fid_msh"]])

    if gmsh_out != 0:
        raise RuntimeError('Gmsh failed in generating this geometry')
    
    fid_dict["fid_xml"] = os.path.join(meshdir, meshfile)
    subprocess.check_output(["dolfin-convert", fid_dict["fid_msh"], fid_dict["fid_xml"]])
    # for debugging:
    # convert2xml(fid_dict["fid_msh"], fid_dict["fid_xml"])
    return dolfin.Mesh(fid_dict["fid_xml"])
    
 
