# (c) Gregor Mitscha-Baude 2016
"2D/3D boxes plus balls"

"""Interface:
A = Box([-1, -1], [1, 1])
m = [0, 0]; r = 0.25
B = Ball(m, r)
C = A | B
C.create_geometry(lc=0.1)
"""
import nanopores.py4gmsh as py4gmsh
from . import box
from .box import (BoxCollection, Float, csgExpression, FacetLoop, Entity,
                BoundaryCollection, union, set_tol)
__all__ = ["Box", "Ball", "EmptySet"]

class BallCollection(BoxCollection):
    def __init__(self, boxes, balls):
        self.balls = list(balls)
        if self.balls:
            self.dim = max(ball.dim for ball in balls)
            self.indexsets = [set() for k in range(self.dim+1)]
        BoxCollection.__init__(self, *boxes)

    def entities_to_gmsh(self, lc=.5, merge=True):
        """
        at this point, all the box-related entities exist and
        all subdomains know the indices of their box boundaries
        DO:
         - compute ball facets and save their indices
         - compute all ball volumes that are not explicitly un-contained in
           domain and whose center lies in domain.
           save them in list / attach to balls.
         - distribute negative indices to subdomains that contain a ball center
           (they can now create their volumes as always)
         - for all subdomains that contain a ball, add their existing surfaces
           as one entry and the faces of their balls as further entries to
           list in self.gmsh_subs.
        """
        # initialize
        self.gmsh_entities = [[None for e in k] for k in self.entities]
        dim = self.dim
        gfacets = self.gmsh_entities[dim-1]

        # compute ball facets and save their indices
        for ball in self.balls:
            #print
            #print ball
            # determine subdomain where ball lies (or none)
            j = ball.boxes[0].indexsets[0].pop() # index of ball center
            sub0 = None
            ball.isubs = []
            ball.ibousubs = []
            for k, sub in enumerate(self.subdomains):
                if j in sub.indexsets[0]:
                    #print "inside", sub.name
                    sub0 = sub

                if (j in sub.indexsets[0] and not sub.csg.excludes(ball) and\
                    not (hasattr(ball, "_added") and ball._added))\
                    or sub.csg.contains(ball):
                    ball.isubs.append(k)
            #if sub0 is None: print "no subdomain"
            # remember to add indices to abstract boundaries whose boundarysub
            # mentions ball
            for k, sub in enumerate(self.boundarysubs):
                if sub.csg.mentions(ball):
                    ball.ibousubs.append(k)
            # remember if ball lies in domain
            if (j in self.indexsets[0] and not self.csg.excludes(ball))\
                or self.csg.contains(ball):
                ball.indomain = True
                #print "in domain"
            else:
                ball.indomain = False
                #print "not in domain"

            # build ball surface in gmsh
            surfs, n = gmsh_ball_surfs(ball, lc)
            # add facets at the end of gmsh_entities[d-1]
            gfacets += surfs # works for lists!
            indices = list(range(len(gfacets)-n, len(gfacets)))
            # add indices and orientations to sub0.bdry() and ball.bdry()
            if sub0 is not None:
                sub0.bdry().indexset |= set(indices)
                for i in indices:
                    sub0.bdry().orients[i] = -1
            # remember indices for ball
            ball.ibdry = indices

        # gmsh all box volumes with balls cut out
        self.entities_to_gmsh_merge(lc)
        gsubs = self.gmsh_subs

        # gmsh all ball volumes that are in domain
        # and add to gmsh_subs where appropriate
        for ball in self.balls:
            if ball.indomain:
                #print "gmshing", ball
                # gmsh ball volume
                subfacets = [gfacets[i] for i in ball.ibdry]
                loop = FacetLoop[dim-1](subfacets)
                gmsh_e = Entity[dim](loop)
                for k in ball.isubs:
                    if gsubs[k] is None:
                        gsubs[k] = gmsh_e
                    elif isinstance(gsubs[k], list):
                        gsubs[k].append(gmsh_e)
                    else:
                        gsubs[k] = [gsubs[k], gmsh_e]
                # update bdry indexsets
                for k in ball.ibousubs:
                    self.boundarysubs[k].bdry().indexsets[dim-1] |= set(ball.ibdry)

        # rebuild boundaries involving balls
        for bou in self.boundaries:
            bou.indexset = bou.csg.evalsets()[dim-1]

    def addsubdomain(self, sub, name):
        assert isinstance(sub, BallCollection)
        sub.name = name
        self.subdomains.append(sub)
        self.boxes = list(set(self.boxes + sub.boxes))
        self.balls = list(set(self.balls + sub.balls))

    def addboundary(self, sub, name):
        assert isinstance(sub, BoxCollection)
        sub.name = name
        self.boundaries.append(sub)
        self.boxes = list(set(self.boxes + sub.boxes))
        # remember all BoundaryCollections involved
        # because we have to compute their indexsets later
        for A in sub.csg.singletons():
            if isinstance(A, BoundaryCollection):
                assert isinstance(A.coll, BallCollection)
                self.boundarysubs.append(A.coll)
                self.boxes = list(set(self.boxes + A.coll.boxes))
                self.balls = list(set(self.balls + A.coll.balls))

    def addball(self, ball, subname="ball", boundaryname="ballb"):
        ball._added = True
        self.addsubdomain(ball, subname)
        self.addboundary(ball.bdry(), boundaryname)

    def addballs(self, balls, subname="ball", boundaryname="ballb"):
        for ball in balls:
            ball._added = True
        sub = union(balls)
        self.addsubdomain(sub, subname)
        self.addboundary(sub.bdry(), boundaryname)

    def _join(self, other):
        boxes = list(set(self.boxes + other.boxes))
        balls = (list(set(self.balls + other.balls))
            if hasattr(other, "balls") else self.balls)
        return BallCollection(boxes, balls)

    def __or__(self, other):
        coll = self._join(other)
        coll.csg = self.csg | other.csg
        return coll

    def __and__(self, other):
        coll = self._join(other)
        coll.csg = self.csg & other.csg
        return coll

    def __sub__(self, other):
        coll = self._join(other)
        coll.csg = self.csg - other.csg
        return coll

class EmptySet(BallCollection):
    def __init__(self, dim=3):
        self.csg = csgExpression(self)
        self.dim = dim
        self.indexsets = [set() for k in range(dim+1)]
        BallCollection.__init__(self, [], [])

    def __repr__(self):
        return "EmptySet()"

class Ball(BallCollection):
    def __init__(self, m, r, lc=None):
        self.m = tuple(Float(x) for x in m)
        self.r = Float(r)
        self.lc = lc
        self.dim = len(m)
        self.csg = csgExpression(self)
        box = Box(m, m)
        BallCollection.__init__(self, [box], [self])

    def __repr__(self):
        return "Ball(%s, %s, %s)" % (self.m, self.r, self.lc)

class Box(BallCollection, box.Box):
    def __init__(self, *args, **params):
        self.balls = []
        box.Box.__init__(self, *args, **params)
        self.indexsets = [set() for k in range(self.dim+1)]


def gmsh_ball_surfs(ball, lc):
    if ball.lc is not None:
        lc=ball.lc
    surfs = add_ball(ball.m, ball.r, lc)
    n = len(surfs)
    return surfs, n

def add_ball(m, r, lc):
    # add ball in 1D, 2D, 3D
    if len(m)==3:
        return py4gmsh.add_ball(m, r, lc, with_volume=False)[2]
    elif len(m)==2:
        return add_circle(m, r, lc)
    elif len(m)==1:
        return [py4gmsh.Point([m[0]-r,0,0]), py4gmsh.Point([m[0]+r,0,0])]
    else:
        raise Exception("Ball center m must have dimension 1, 2 or 3.")

def add_circle(m, r, lc):
    m0, m1 = m[0], m[1]
    point = lambda x, y: py4gmsh.Point([x,y,0], lc)
    # add points.
    p = [point(m0, m1), point(m0+r, m1), point(m0, m1+r),
         point(m0-r, m1), point(m0, m1-r)]
    # add circle lines
    return [py4gmsh.Circle([p[1], p[0], p[2]]),
            py4gmsh.Circle([p[2], p[0], p[3]]),
            py4gmsh.Circle([p[3], p[0], p[4]]),
            py4gmsh.Circle([p[4], p[0], p[1]])]

if __name__ == "__main__":
    # unit square
    A = Box([-1, -1], [1, 1]) | Box([0.5, -1.2], [1.2, -0.5])
    B = Ball([2, 0], 0.5, lc=1.)
    balls = [Ball((-0.5, -0.5), 0.4), Ball((0.5, 0.5), 0.4)]
    smallballs = [Ball((x*.1, -x*.1), r=0.08, lc=0.02) for x in range(-8,12,2)]
    # union
    C = A | B
    C.addsubdomains(box=A-B, circle=B)
    C.addboundaries(boxb=(A | B).boundary())
    C.addballs(balls) # adds ball boundaries automatically
    C.addballs(smallballs, "small", "smallb")

    C.create_geometry(lc=0.1)
    print(C.geo)
    C.plot()