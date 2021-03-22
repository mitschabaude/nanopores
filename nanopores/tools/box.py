'''
python class to handle n-dimensional rectangular boxes and unions of them

note: very naive implementation of union, suitable for a couple 100 boxes
'''

# TODO: implement tests for points inside (sub-)domain by unpacking synonymes
#       and traversing expression tree.
# TODO: customizable length scales
# TODO: to make merge always applicable, would need check for connectedness
#       and handle Surface in Volume cases
# TODO: it would be nice to have some sort of "interactive mode" where you add
#       domains and boundaries in order and they are plotted immediately
#       -- for showing off.
#       this could use the module.__getattr__ trick (stackoverflow)
#       or be implemented properly with update of entities etc.
# TODO: reasonable verbosity

import sys
import os
from itertools import product, combinations
import nanopores.py4gmsh as py4gmsh
import dolfin
from functools import reduce

def printnow(s):
    print(s, end=' ')
    sys.stdout.flush()

__all__ = ["BoxCollection", "Box", "Interval", "Log"]
MESHDIR = "/tmp/nanopores"

class BoxCollection(object):
    """
    collection of boxes together with constructive solid geometry (CSG)
    expression tree."""

    def __init__(self, *boxes):
        self.boxes = list(boxes)
        self.subdomains = []
        self.boundaries = []
        self.boundarysubs = []
        self.indexset = set()
        self.params = {}
        #self.csg = union(map(csgExpression, boxes))

    def compute_entities(self):
        #print self.boxes
        dic = multi_box_union(self.boxes)
        self.__dict__.update(dic)
        #self.indexset = self._csg().eval()
        self.indexsets = self.csg.evalsets()
        d = self.dim
        for sub in self.subdomains + self.boundarysubs:
            sub.indexsets = sub.csg.evalsets()
            sub.indexset = sub.indexsets[d] & self.indexsets[d]
        # make sure that subdomains cover all of domain:
        if self.subdomains:
            rest = self - union(self.subdomains)
            rest.indexsets = rest.csg.evalsets()
            rest.indexset = rest.indexsets[d]
            if rest.indexset:
                self.addsubdomain(rest, "rest")
        else:
            # create dummy subdomain
            self.addsubdomain(self, "domain")
            self.indexset = self.indexsets[d]
        if not self.boundaries:
            # create dummy domain boundary
            self.addboundary(Boundary(self), "boundary")
            self.indexset = self.indexsets[d]
        #print "BOUNDARIES", [b.name for b in self.boundaries]
        #print self.boundarysubs
        #print self.boundary()

    def compute_boundaries(self, merge=True):
        # call after compute_entities
        # need list of domains of which to compute boundary
        d = self.dim

        subs = self.boundarysubs
        if merge:
            subs += self.subdomains
        orients = [eval(o+"1") for o in _orientations(d-1)]
        #print orients
        for sub in subs:
            #print sub.name, ":"
            iset = set()
            odict = dict()
            for i in sub.indexset:
                e = self.entities[d][i]
                for f, o in zip(_facets(e), orients):
                    iface = self.entities[d-1].index(f)
                    #print odict
                    #print iset
                    #print iface
                    if iface in iset:
                        iset.remove(iface)
                        odict[iface] = odict[iface] + o
                    else:
                        iset.add(iface)
                        odict[iface] = o
            sub.bdry().indexset = iset
            sub.bdry().indexsets = [set() for i in range(d+1)]
            sub.bdry().indexsets[d-1] = iset
            sub.bdry().orients = odict

        for sub in self.boundaries:
            # ATTENTION: if boundary is not contained in domain, its entities
            # will still be created. this is addressed by deleting Nones in the
            # gmsh surfaces defining the physical boundary, which can silently
            # cause an empty boundary in Geometry and e.g. BCs without effect.
            #sub.indexset = sub.csg.eval() #& self.indexsets[d-1] # <-- wrong
            sub.indexsets = sub.csg.evalsets()
            sub.indexset = sub.indexsets[d-1]

    def entity_to_gmsh(self, e, dim, lc, gmshself=True):
        # do not duplicate entity in gmsh
        i = self.entities[dim].index(e)
        gmsh_e = self.gmsh_entities[dim][i]
        if gmsh_e is not None:
            return gmsh_e

        if dim==0: # create Point
            e = e + tuple(0. for i in range(3 - self.dim))
            gmsh_e = py4gmsh.Point(e, lc)
            self.gmsh_entities[0][i] = gmsh_e
            #print gmsh_e, e
            return gmsh_e

        # dim>0: recursively generate facets and entity itself
        facets = _facets(e)
        facets = [self.entity_to_gmsh(f, dim-1, lc)
            for f in facets]
        orient = _orientations(dim-1)
        loop = FacetLoop[dim-1]([o+s for o, s in zip(orient, facets)])
        if gmshself:
            gmsh_e = Entity[dim](loop)
            self.gmsh_entities[dim][i] = gmsh_e
            #print gmsh_e, e
            return gmsh_e

    def entities_to_gmsh_nomerge(self, lc):
        dim = self.dim
        # add facets of entities of full dimension
        for i, e in enumerate(self.entities[dim]):
            if not i in self.indexsets[dim]:
                continue
            self.entity_to_gmsh(e, dim, lc)

    def entities_to_gmsh_merge(self, lc):
        dim = self.dim
        # build full subdomains and their boundaries
        facets = self.entities[dim-1]
        gmsh = self.gmsh_entities[dim-1]
        self.gmsh_subs = [None for sub in self.subdomains]
        for j, sub in enumerate(self.subdomains):
            I = sub.bdry().indexset
            # gmsh won't compile empty volume
            if not I: continue
            # gmsh facets
            subfacets = []
            for i in I:
                # do not duplicate entity in gmsh
                if gmsh[i] is None:
                    self.entity_to_gmsh(facets[i], dim-1, lc)
            # gmsh sub
            orients = sub.bdry().orients
            dic = {1:"+", -1:"-"}
            subfacets = [dic[orients[i]] + gmsh[i] for i in I]
            loop = FacetLoop[dim-1](subfacets)
            gmsh_e = Entity[dim](loop)
            self.gmsh_subs[j] = gmsh_e

    def entities_to_gmsh(self, lc=.5, merge=True):
        # initialize
        self.gmsh_entities = [[None for e in k] for k in self.entities]
        if merge:
            self.entities_to_gmsh_merge(lc)
        if not merge:
            self.entities_to_gmsh_nomerge(lc)

    def physical_to_gmsh(self, merge=True):
        # call after entities_to_gmsh
        self.dimt = dimt = -1 + sum(1 for en in self.entities if en)

        if merge:
            for sub, vol in zip(self.subdomains, self.gmsh_subs):
                if vol is not None:
                    py4gmsh.PhysicalVolume(vol, sub.name, dimt)
                else:
                    py4gmsh.NoPhysicalVolume(sub.name)

        else:
            for sub in self.subdomains:
                vols = [self.gmsh_entities[dimt][i] for i in sub.indexset]
                py4gmsh.PhysicalVolume(vols, sub.name, dimt)

        for sub in self.boundaries:
            surfs = [self.gmsh_entities[dimt-1][i] for i in sub.indexset]
            surfs = [s for s in surfs if s is not None]
            if len(surfs)>0:
                py4gmsh.PhysicalSurface(surfs, sub.name, dimt)
            else:
                py4gmsh.NoPhysicalSurface(sub.name)

    def create_geometry(self, lc=.5, merge=True):
        with Log("computing geometrical entities..."):
            self.compute_entities()

        with Log("computing boundaries..."):
            self.compute_boundaries(merge)

        with Log("computing gmsh entities..."):
            self.entities_to_gmsh(lc, merge)
            self.physical_to_gmsh(merge)
        #entities_to_gmsh(self.entities, self.indexsets, self.esets, lc=lc)
        #self.gmsh_entities = gmsh_entities

        self.geo = to_mesh()
        self.geo.params = self.params
        if hasattr(self, "synonymes"):
            self.geo.import_synonymes(self.synonymes)
        return self.geo

    def write_gmsh_code(self, lc=1., merge=True):
        with Log("computing geometrical entities..."):
            self.compute_entities()

        with Log("computing boundaries..."):
            self.compute_boundaries(merge)

        with Log("computing gmsh entities..."):
            self.entities_to_gmsh(lc, merge)
            self.physical_to_gmsh(merge)

    def insert_points(self, points, lc=1., dim=None, forbidden=()):
        if not points: return
        if dim is None: dim = len(points[0])
        vol = {1: "Line", 2: "Surface", 3: "Volume"}[dim]
        with Log("inserting points..."):
            for x in points:
                x_ = [x[i] if i<dim else 0. for i in range(3)]
                surf, name = self.get_gmsh_sub(x)
                #print p, surf, name
                if name not in forbidden:
                    p = py4gmsh.Point(x_, lc)
                    py4gmsh.raw_code(["Point{%s} In %s{%s};" %(p, vol, surf)])

    def get_gmsh_sub(self, x):
        d = self.dim
        for j, sub in enumerate(self.subdomains):
            for i in sub.indexset:
                e = self.entities[d][i]
                if _inside_entity(x, e):
                    #print "x in", sub.name
                    return self.gmsh_subs[j], sub.name
        return None, None

    def code_to_mesh(self):
        self.geo = to_mesh()
        self.geo.params = self.params
        if hasattr(self, "synonymes"):
            self.geo.import_synonymes(self.synonymes)
        return self.geo

    def recreate_geometry(self):
        self.geo = geo_from_meshdir()
        self.geo.params = self.params
        if hasattr(self, "synonymes"):
            self.geo.import_synonymes(self.synonymes)
        return self.geo

    def plot(self, sub=False):
        geo = self.create_geometry() if not hasattr(self, "geo") else self.geo
        if hasattr(geo, "subdomains"):
            dolfin.plot(geo.subdomains)
            dolfin.plot(geo.boundaries)
        else:
            dolfin.plot(geo.mesh)
        dolfin.interactive()

    def addsubdomain(self, sub, name):
        assert isinstance(sub, BoxCollection) or isinstance(sub, Box)
        sub.name = name
        self.subdomains.append(sub)
        self.boxes = list(set(self.boxes + sub.boxes))
        #self.boxes = _unique(self.boxes + sub._boxes())

    def getsubdomain(self, name):
        # helper function to find subdomain from its name
        for sub in self.subdomains:
            if sub.name == name:
                return sub

    def addsubdomains(self, **subdomains):
        for name, sub in list(subdomains.items()):
            self.addsubdomain(sub, name)

    def addboundary(self, sub, name):
        assert isinstance(sub, BoxCollection) or isinstance(sub, Box)
        sub.name = name
        self.boundaries.append(sub)
        self.boxes = list(set(self.boxes + sub.boxes))
        # remember all BoundaryCollections involved
        # because we have to compute their indexsets later
        for A in sub.csg.singletons():
            if isinstance(A, BoundaryCollection):
                self.boundarysubs.append(A.coll)
                self.boxes = list(set(self.boxes + A.coll.boxes))

    def getboundary(self, name):
        for sub in self.boundaries:
            if sub.name == name:
                return sub

    def addboundaries(self, **boundaries):
        for name, sub in list(boundaries.items()):
            self.addboundary(sub, name)

    def bdry(self):
        return Boundary(self)

    def boundary(self):
        return Boundary(self)

    def copy(self):
        # copy constructor
        newcol = eval(repr(self))
        # reconstruct subdomains
        for sub in self.subdomains:
            newsub = eval(repr(sub))
            newcol.addsubdomain(newsub, sub.name)
        for sub in self.boundaries:
            newsub = eval(repr(sub))
            newcol.addboundary(newsub, sub.name)
        # copy other stuff
        newcol.params = dict(self.params)
        #newcol.indexset = set(self.indexset)
        if hasattr(self, "synonymes"):
            newcol.synonymes = dict(self.synonymes)
        return newcol

    # TODO bad hack??
    # list of class names that are allowed to override set operators
    operator_priority = ["BallCollection", "Ball"]

    def __or__(self, other):
        if type(other).__name__ in self.operator_priority:
            return other | self
        boxes = list(set(self.boxes + other.boxes))
        coll = BoxCollection(*boxes)
        coll.csg = self.csg | other.csg
        return coll

    def __and__(self, other):
        if type(other).__name__ in self.operator_priority:
            return other & self
        boxes = list(set(self.boxes + other.boxes))
        coll = BoxCollection(*boxes)
        coll.csg = self.csg & other.csg
        return coll

    def __sub__(self, other):
        if type(other).__name__ in self.operator_priority:
            return other - self
        boxes = list(set(self.boxes + other.boxes))
        coll = BoxCollection(*boxes)
        coll.csg = self.csg - other.csg
        return coll

    def __contains__(self, coll):
        return self.csg.contains(coll)

    #def __str__(self):
    #    return str(self.csg)

    def __repr__(self):
        return str(self.csg)

def Interval(a, b):
    # simple wrapper for Box
    return Box([a], [b])

class Box(BoxCollection):

    def __init__(self, a=None, b=None, intervals=None, center=None,
                 l=None, w=None, h=None):
        if intervals is not None:
            a = tuple(Float(x) for x,y in intervals)
            b = tuple(Float(y) for x,y in intervals)
        if center is not None:
            lengths = [l, w, h][:len(center)]
            a = [c - x*.5 for c,x in zip(center, lengths)]
            b = [c + x*.5 for c,x in zip(center, lengths)]
        # a, b should be tuples
        assert len(a) == len(b)
        a, b = _sort(a, b)
        a = tuple(Float(x) for x in a)
        b = tuple(Float(x) for x in b)
        self.intervals = list(zip(a, b))
        self.dim = len(a)
        self.a = a
        self.b = b

        # determine topological dimension
        self.dimt = 0
        for i in self.intervals:
            if not i[0] == i[1]:
                self.dimt += 1

        self.csg = csgExpression(self)
        BoxCollection.__init__(self, self)

    #def __str__(self):
    #    return "Box(%s)" %("x".join(["[%s, %s]" %i for i in self.intervals]),)

    def __repr__(self):
        return "Box(%s, %s)" %(self.a, self.b)

    def __cmp__(self, other):
        return cmp(self.intervals, other.intervals)

    def __iter__(self):
        return iter(self.intervals)

    def __getitem__(self, i):
        return self.intervals[i]

    def boundary(self, *strings):
        facets = _facets(self.intervals)
        if strings:
            ind = [_boundaries[self.dim][s] for s in strings]
            facets = [facets[i] for i in ind]
        boxes = []
        for facet in facets:
            intervals = [(f if isinstance(f, tuple) else (f,f)) for f in facet]
            boxes.append(Box(intervals=intervals))
        return union(boxes)

class BoundaryCollection(BoxCollection):
    "boundary of a BoxCollection, which is a BoxCollection itself"
    """this is just an empty collection with csg initialized as singleton,
    which means self.indexset has to be built manually."""

    def __init__(self, coll):
        self.csg = csgExpression(self)
        self.coll = coll
        if hasattr(coll, "dim"):
            self.indexsets = [set() for i in range(coll.dim+1)]
        BoxCollection.__init__(self)

    def __repr__(self):
        return "Boundary(%s)" % (repr(self.coll),)
    def __str__(self):
        return "Boundary(%s)" % (repr(self.coll),)

def Boundary(coll):
    """wrapper for BoundaryCollection that refers to a single object, i.e.
    Boundary(A) == Boundary(A)"""
    if not hasattr(coll, "_bdry"):
        coll._bdry = BoundaryCollection(coll)
    return coll._bdry


_boundaries = {
1 : dict(
    left = 0,
    right = 1,
),
2 : dict(
    left = 0,
    right = 1,
    top = 3,
    bottom = 2,
),
3 : dict(
    left = 0,
    right = 1,
    front = 2,
    back = 3,
    top = 5,
    bottom = 4,
)}

class csgExpression(object):
    # either singleton (box) or made up of two csgExpressions via set operations

    def __init__(self, *args):
        if len(args) == 1:
            self.sginit(*args)
        else:
            self.opinit(*args)

    def sginit(self, A):
        self.A = A
        self.singleton = True

    def opinit(self, op, A, B):
        self.op = op
        self.A = A
        self.B = B
        self.singleton = False

    def __repr__(self):
        if self.singleton:
            return repr(self.A)
        else:
            return "(%s %s %s)" %(repr(self.A), self.op, repr(self.B))

    def eval(self): # is designed to return a set
        if self.singleton:
            return self.A.indexset
        if self.op == "|":
            return self.A.eval() | self.B.eval()
        elif self.op == "&":
            return self.A.eval() & self.B.eval()
        elif self.op == "-":
            return self.A.eval() - self.B.eval()

    # TODO: get rid of eval(). having an indexset for every dimension is better
    def evalsets(self): # indexset for every dimension
        if self.singleton:
            return self.A.indexsets
        if self.op == "|":
            return [a | b for a, b in zip(self.A.evalsets(), self.B.evalsets())]
        elif self.op == "&":
            return [a & b for a, b in zip(self.A.evalsets(), self.B.evalsets())]
        elif self.op == "-":
            return [a - b for a, b in zip(self.A.evalsets(), self.B.evalsets())]

    def singletons(self):
        "traverse expression tree and return list of leafs"
        if self.singleton:
            return [self.A]
        else:
            return self.A.singletons() + self.B.singletons()

    def contains(self, C):
        # this is only if C is singleton and tests containment as if all leafs
        # were disjoint. useful for entities like balls.
        if self.singleton:
            return True if self.A is C else False
        elif self.op == "|":
            return self.A.contains(C) or self.B.contains(C)
        elif self.op == "&":
            return self.A.contains(C) and self.B.contains(C)
        elif self.op == "-":
            return self.A.contains(C) and not self.B.contains(C)

    def mentions(self, C):
        return any(C is A for A in self.singletons())

    def excludes(self, C):
        return self.mentions(C) and not self.contains(C)

    # boolean operations
    def __or__(self, other):
        return csgExpression("|", self, other)
    def __and__(self, other):
        return csgExpression("&", self, other)
    def __sub__(self, other):
        return csgExpression("-", self, other)

# TOL = None # much faster
TOL = 1e-10 # more secure

def set_tol(tol):
    global TOL
    TOL = tol

def Float(a):
    if TOL is None:
        return a
    else:
        return ExactFloat(a)

class ExactFloat(float):

    def __cmp__(self, other):
        return _cmp(self, other, TOL)

    # TODO: this is performance-critical, so minimize function calls
    def __eq__(self, other):
        if isinstance(other, float):
            return abs(self - other) < TOL
            #return cmp(self, other) == 0
        else:
            return False

    def __ne__(self, other):
        return not self == other

def union(seq):
    try:
        return reduce(lambda x, y: x | y, seq)
    except TypeError:
        return set()

def intersection(seq):
    return reduce(lambda x, y: x & y, seq)

def _sort(a, b):
    a_ = tuple(map(min, zip(a,b)))
    b_ = tuple(map(max, zip(a,b)))
    return a_, b_

def _cmp(s, t, tol):
    " test floats for approximate equality "
    if s > t+tol:
        return 1
    elif s < t-tol:
        return -1
    else:
        return 0

def _unique(seq):
    # version of list(set(seq)) but using list membership test
    # (weaker than set membership test, for approximate Float equality)
    unique = []
    [unique.append(x) for x in seq if not unique.count(x)]
    return unique

def _inside_entity(x, entity):
    "test whether point x = (x0,x1,x2) lies in entity"
    for t, e in zip(x, entity):
        if isinstance(e, tuple):
            if not e[0] <= t <= e[1]:
                return False
        else:
            if not e == t:
                return False
    return True

def multi_interval_union(intvs):
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
    subs = list(zip(points[:-1], points[1:]))
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


def _map(f, seq):
    """ version of map() for function with multiple return values.
    instead of a list of tuples (like map), return a tuple of lists. """
    return tuple(map(list, list(zip(*list(map(f, seq))))))

def multi_box_union(boxes): #, facets=[]):
    " return union of boxes as collection of disjoint boxes "
    # TODO: atm this modifies the boxes with an "indexset" attribute
    #       which doesn't seem perfectly clean
    allboxes = boxes #+ facets

    # get dimension; all boxes must have the same
    dim = boxes[0].dim
    assert all([dim == box.dim for box in boxes])

    # get list of disjoint intervals for every dimension
    nodes, nsets, intvs, isets = _map(multi_interval_union, zip(*allboxes))

    D = list(range(dim)) # [0,...,dim-1]
    D1 = list(range(dim+1)) # [0,..,dim]
    entities = [[] for k in D1] # entity := d-tuple of tuples/floats <=> box
    esets = [[] for k in D1] # eset := {indices of the parent boxes of an entity}
    for box in allboxes:
        box.indexsets = [set() for i in D1]
        #box.indexset = set()

    # collect all admissible k-entities for k = 1,...,d
    # 0-entity = vertex, 1-entity = edge, ...
    for k in D1:
        j = 0
        for tup in combinations(D, k):
            k_entities = product(*[(intvs[i] if i in tup else nodes[i]) for i in D])
            k_esets1D = product(*[(isets[i] if i in tup else nsets[i]) for i in D])
            for entity, eset1D in zip(k_entities, k_esets1D):
                eset = set.intersection(*eset1D)
                if eset:
                    entities[k].append(entity)
                    esets[k].append(eset)
                    for ii in eset:
                        allboxes[ii].indexsets[k].add(j)
                    j += 1 # counts all k-entities

        #print "k=%s:"%k, entities[k]
        #print
    # choose correct indexset for every box
    for box in allboxes:
        box.indexset = box.indexsets[box.dimt]
        #print box
        #print box.indexsets

    return dict(nodes=nodes, entities=entities, esets=esets, dim=dim)

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

# DEPRECATED
def entities_to_gmsh(entities, indexsets, esets, lc=.5):

    global gmsh_entities
    global dim
    global dimt

    dim = len(entities)-1
    dimt = -1
    for en in entities:
        if en:
            dimt += 1

    gmsh_entities = [[None for e in k] for k in entities]

    # a shortcut
    _gmsh = lambda entity, k: gmsh_entities[k][entities[k].index(entity)]

    # add points
    for i, e in enumerate(entities[0]):
        #print e
        #print esets[0][i]
        # py4gmsh.Point expects three values x,y,z
        e = e + tuple(0. for i in range(3 - dim))
        # create point
        gmsh_entities[0][i] = py4gmsh.Point(e, lc)

    # add entities of dimension > 0
    for k in range(1, dim+1):
        orient = _orientations(k-1)

        for i, e in enumerate(entities[k]):
            #print e
            #print esets[k][i]

            # skip code creation for entity if not part of domain
            # what do we do with unnessecary facets??
            if k==dim and i not in indexsets[k]:
                #print "NOT INSIDE"
                #print
                continue

            # get surrounding facets (x,y,z in this order)
            facets = _facets(e)

            #for f in facets:
                #print f,

            # find corresponding gmsh facets
            facets = [_gmsh(l, k-1) for l in facets]

            # create facet loop with correct orientations
            #print facets
            loop = FacetLoop[k-1]([o+s for o, s in zip(orient, facets)])
            #print

            # add entity
            gmsh_entities[k][i] = Entity[k](loop)

    return gmsh_entities

def to_mesh(clscale=1., pid=""):
    pid = str(os.getpid())
    with Log("executing gmsh..."):
        py4gmsh.raw_code(["General.ExpertMode = 1;"])
        py4gmsh.raw_code(["Mesh.Algorithm3D = 2;"])
        code = py4gmsh.get_code()
        meta = py4gmsh.get_meta()
        if not meta["physical_domain"]:
            pass

        import subprocess
        import nanopores
        inputfile = "input%s.geo" %pid
        outfile = "out%s.msh" %pid
        meshfile = "mesh%s.xml" %pid

        # create path/to/nanoporesdata/gid/mesh if not already there
        meshdir = MESHDIR
        if not os.path.exists(meshdir):
            os.makedirs(meshdir)

        xml_sub = meshdir+"/mesh%s_physical_region.xml" %pid
        xml_bou = meshdir+"/mesh%s_facet_region.xml" %pid
        if os.path.exists(xml_sub): os.remove(xml_sub)
        if os.path.exists(xml_bou): os.remove(xml_bou)

        fid_dict = {"fid_geo": os.path.join(meshdir, inputfile),
                    "fid_msh": os.path.join(meshdir, outfile)}

        # save code to .geo file
        fobj = open(fid_dict["fid_geo"], "w")
        fobj.write(code)
        fobj.close()

        # after writing the geo file, call gmsh
        gmsh_out = subprocess.call(["gmsh", "-3", "-v", "1","-clscale", "%f" %clscale,
                         fid_dict["fid_geo"], "-o", fid_dict["fid_msh"], "-optimize"])

        if gmsh_out != 0:
            raise RuntimeError('Gmsh failed in generating this geometry')

    with Log("converting to dolfin..."):
        fid_dict["fid_xml"] = os.path.join(meshdir, meshfile)
        subprocess.check_output(["dolfin-convert", fid_dict["fid_msh"], fid_dict["fid_xml"]])
        # for debugging:
        # convert2xml(fid_dict["fid_msh"], fid_dict["fid_xml"])
        mesh = dolfin.Mesh(fid_dict["fid_xml"])

    #print(meta)
    with open('%s/%s.txt' % (meshdir, "meta%s" %pid), 'w') as f:
        f.write(repr(meta))

    physdom = meta.pop("physical_domain")
    physbou = meta.pop("physical_boundary")
    subdomains = dolfin.MeshFunction("size_t", mesh, xml_sub) if physdom else None
    boundaries = dolfin.MeshFunction("size_t", mesh, xml_bou) if physbou else None

    return nanopores.Geometry(None, mesh, subdomains, boundaries, physdom, physbou)


def geo_from_meshdir(DIR=MESHDIR):
    import nanopores
    mesh = dolfin.Mesh(DIR+"/mesh.xml")
    subdomains = dolfin.MeshFunction("size_t", mesh, DIR+"/mesh_physical_region.xml")
    boundaries = dolfin.MeshFunction("size_t", mesh, DIR+"/mesh_facet_region.xml")

    with open(DIR+"/meta.txt", "r") as f:
        meta = eval(f.read())

    physdom = meta.pop("physical_domain")
    physbou = meta.pop("physical_boundary")

    return nanopores.Geometry(None, mesh, subdomains, boundaries, physdom, physbou)


class Log(object):
    def __init__(self, msg):
        self.msg = msg
    def __enter__(self):
        printnow(self.msg)
        dolfin.tic()
    def __exit__(self, *args):
        print("%.2g s" %(dolfin.toc(),))