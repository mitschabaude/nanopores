"2D/3D boxes plus balls"

"""Interface:
A = Box([0, 0], [1, 1])
m = [0, 0]; r = 0.25
B = Ball(m, r)
C = A | B
C.create_geometry(lc=0.1)
"""
import nanopores.py4gmsh as py4gmsh
from box import BoxCollection, Box, Float, csgExpression

# TODO gmsh_ball_surfs for 1D, 2D
# TODO make ball union domain possible -> gmsh_sub with multiple vols
# TODO no logic for balls yet, i.e. balls get cut out of domains they are
#      contained in.

class BallCollection(BoxCollection):
    def __init__(self, boxes, balls):
        self.balls = list(balls)
        if self.balls:
            self.dim = max(ball.dim for ball in balls)
            self.indexsets = [set() for k in range(self.dim+1)] 
        BoxCollection.__init__(self, *boxes)
    
    def entities_to_gmsh(self, lc=.5, merge=True):        
        # initialize
        self.gmsh_entities = [[None for e in k] for k in self.entities]
        gfacets = self.gmsh_entities[self.dim-1]
        
        for ball in self.balls:
            # determine subdomain where ball lies (or none)
            j = ball.boxes[0].indexsets[0].pop() # index of ball center
            sub0 = None
            for sub in self.subdomains:
                if j in sub.indexsets[0]:
                    print sub.name
                    sub0 = sub
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
            ball.bdry().indexset |= set(indices)
            if not hasattr(ball.bdry(), "orients"):
                ball.bdry().orients = dict()
            for i in indices:
                ball.bdry().orients[i] = 1
        
        # TODO is this necessary??
        for sub in self.subdomains:
            for singleton in sub.csg.singletons():
                if isinstance(singleton, Ball):
                    print sub.name
                    iset = singleton.bdry().indexset
                    sub.bdry().indexset |= iset
                    for i in iset:
                        sub.bdry().orients[i] = 1
                    
        # gmsh everything else
        self.entities_to_gmsh_merge(lc)
        # rebuild boundaries involving balls
        for bou in self.boundaries:
            bou.indexset = bou.csg.evalsets()[self.dim-1]
        
    def _join(self, other):
        boxes = list(set(self.boxes + other.boxes))
        balls = list(set(self.balls + other.balls)) if hasattr(other, "balls")\
            else self.balls
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
        
def gmsh_ball_surfs(ball, lc):
    if ball.lc is not None:
        lc=ball.lc
    surfs = py4gmsh.add_ball(ball.m, ball.r, lc, with_volume=False)[2]
    n = len(surfs)
    return surfs, n
    
if __name__ == "__main__":
    # unit square
    A = Box([-1, -1, -1], [1, 1, 1])
    # circle with r=0.5 centered at origin
    B = Ball([0, 0, 0], 0.5)
    B1 = Ball([2, 0, 0], 0.5, lc=0.05)
    # union
    C = A | B | B1
    C.addsubdomains(box=A, ball=B, ball1=B1)
    #C.addsubdomains(box=A, balls=(B|B1)) #, ball=B, ball1=B1)
    C.addboundaries(boxb=A.boundary()) #, ballb=B.boundary())#, ball1b=B1.boundary())
    
    C.create_geometry(lc=0.2)
    print C.geo
    from nanopores import plot_sliced
    plot_sliced(C.geo)
    C.plot()