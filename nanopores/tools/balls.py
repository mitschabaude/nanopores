"2D/3D boxes plus balls"

"""Interface:
A = Box([0, 0], [1, 1])
m = [0, 0]; r = 0.25
B = Ball(m, r)
C = A | B
C.create_geometry(lc=0.1)
"""
import nanopores.py4gmsh as py4gmsh
import box
from box import (BoxCollection, Float, csgExpression, FacetLoop, Entity,
                BoundaryCollection, union)

# TODO gmsh_ball_surfs for 1D, 2D
# FIXME ball boundaries do not seem to work

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
            gfacets += surfs
            indices = range(len(gfacets)-n, len(gfacets))
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
        self.addboundary(ball.boundary(), boundaryname)
        
    def addballs(self, balls, subname="ball", boundaryname="ballb"):
        for ball in balls:
            ball._added = True
        sub = union(balls)
        self.addsubdomain(sub, subname)
        self.addboundary(sub.boundary(), boundaryname)
        
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
    surfs = py4gmsh.add_ball(ball.m, ball.r, lc, with_volume=False)[2]
    n = len(surfs)
    return surfs, n
    
if __name__ == "__main__":
    # unit square
    A = Box([-1]*3, [1]*3) | Box([-1.5, -1, -1.5], [-0.95, 1, -0.95])
    B1 = Ball([2, 0, 0], 0.5, lc=0.05)
    B = [Ball((-0.5, 0, -0.5), 0.4), Ball((0.5, 0, 0.5), 0.4)]
    # union
    C = A | B1
    C.addsubdomains(box=A|B1)
    C.addballs(B)
    C.addboundaries(boxb=(A | B1).boundary())
    
    C.create_geometry(lc=0.1)
    print C.geo
    from nanopores import plot_sliced
    from dolfin import interactive, plot
    plot_sliced(C.geo)
    plot(C.geo.boundaries, title="boundaries")
    interactive()