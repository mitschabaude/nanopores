# (c) 2017 Gregor Mitscha-Baude
"""this was intended to perform boolean operations on polygons,
but this seems to be too hard without relying on external software.
now it is simply for preparing the polygons involved in a pore geometry where
only the pore protein crosssection and parameters have to be plugged in."""
import numpy as np
from bisect import bisect
from collections import OrderedDict
from matplotlib import pyplot as plt
import matplotlib.path as mpath
from nanopores.tools.utilities import Params

# convention:
# a - b - c - d numbering of polygon corners
# start from lower left in clockwise direction
# in general, polygon nodes are ordered in CLOCKWISE direction (cccw means hole)
def irange(i0, i1, n, sign=1):
    # index range modulo n
    l = [i0]
    i = i0
    while i != i1:
        i = (i + sign) % n
        l.append(i)
    return l

class Polygon(object):
    TOL = 0.1

    def __init__(self, nodes):
        self.nodes = [(float(x[0]), float(x[1])) for x in nodes]
        self.init_edges()

    def __repr__(self):
        return "Polygon(%s)" % (self.nodes,)

    def __contains__(self, x):
        return winding_number(self, np.array([x])) > 0

    def copy(self):
        return Polygon(self.nodes)

    def init_edges(self):
        nodes = self.nodes
        self.edges = zip(nodes, nodes[1:] + nodes[0:1])

    def inside(self, x, radius=0., newpath=False):
        if newpath or not hasattr(self, "path"):
            self.path = mpath.Path(np.array(self.nodes[::-1]), closed=True)
        #return np.array([self.path.contains_point(t) for t in x])
        return self.path.contains_points(x, radius=radius)

    def inside_single(self, x, radius=0., newpath=False):
        x = np.array([np.sqrt(x[0]**2 + x[1]**2), x[2]])
        if newpath or not hasattr(self, "path"):
            self.path = mpath.Path(np.array(self.nodes[::-1]), closed=True)
        return self.path.contains_point(x, radius=radius)
        #return self.path.contains_points(np.array([x]), radius=radius)[0]
        
    def inside_winding(self, r, z):
        rz = np.column_stack([r, z])
        return winding_number(self, rz) > 0

    def intersections(self, z, axis=1):
        z = float(z)
        return {v: vnodes for edge in self.edges \
                          for v, vnodes in self.intersect_edge(edge, z, axis)}

    def left_intersection(self, z):
        I = self.intersections(z)
        if not I: return
        x = min(I.keys(), key=lambda t: t[0])
        self.add(x, I[x])
        return x

    def right_intersection(self, z):
        I = self.intersections(z)
        if not I: return
        x = max(I.keys(), key=lambda t: t[0])
        self.add(x, I[x])
        return x

    def top_intersection(self, r):
        I = self.intersections(r, 0)
        if not I: return
        x = max(I.keys(), key=lambda t: t[1])
        self.add(x, I[x])
        return x

    def bottom_intersection(self, r):
        I = self.intersections(r, 0)
        if not I: return
        x = min(I.keys(), key=lambda t: t[1])
        self.add(x, I[x])
        return x
    
    def add_intersection(self, side, t):
        # replaces four functions above
        # side in {"top", "bottom", "left", "right"}
        axis = 0 if side in {"top", "bottom"} else 1
        ext = min if side in {"bottom", "left"} else max
        I = self.intersections(t, axis=axis)
        if not I: return
        x = ext(I.keys(), key=lambda t: t[1 - axis])
        self.add(x, I[x])
        return x

    def all_intersections(self, z, axis=1):
        I = self.intersections(z, axis)
        todo = {k: v for k, v in I.items() if len(v)>0}
        while todo:
            x = todo.keys()[0]
            self.add(x, todo[x])
            #print "PROCESSING", x, todo[x]
            I = self.intersections(z, axis)
            todo = {k: v for k, v in I.items() if len(v)>0}

        return I.keys()

    def intersect_edge(self, edge, z, axis=1):
        if axis == 0:
            return self.intersect_edge_vertical(edge, z)
        # return intersection points with line z=z
        # and vertices to provide context for insertion
        x, y = edge[0], edge[-1]
        x0, x1 = x
        y0, y1 = y
        # special case
        if x1 == y1 == z:
            return {(x, ()), (y, ())}
        elif x1 <= z < y1 or y1 < z <= x1:
            # special case
            if abs(y1 - x1) < 1e-5*self.TOL:
                return {(x, ()), (y, ())}
            # x1 + (y1 - x1)*t == z
            # at this point uniquely solvable in (0,1)
            t = (z - x1)/(y1 - x1)
            v = (x0 + (y0 - x0)*t, z)
            if self.close(x, v):
                vnodes = () if x == v else (x,)
            elif self.close(y, v):
                vnodes = () if y == v else (y,)
            else:
                vnodes = (x, y)
            return {(v, vnodes)}
        else:
            return set()

    def intersect_edge_vertical(self, edge, r):
        # return intersection points with line r=r
        # and vertices to provide context for insertion
        x, y = edge
        x0, x1 = x
        y0, y1 = y
        # special case
        if x0 == y0 == r:
            return {(x, ()), (y, ())}
        elif x0 <= r < y0 or y0 < r <= x0:
            # special case
            if abs(y0 - x0) < 1e-5*self.TOL:
                return {(x, ()), (y, ())}
            # x0 + (y0 - x0)*t == r
            # at this point uniquely solvable in (0,1)
            t = (r - x0)/(y0 - x0)
            v = (r, x1 + (y1 - x1)*t)
            if self.close(x, v):
                vnodes = () if x == v else (x,)
            elif self.close(y, v):
                vnodes = () if y == v else (y,)
            else:
                vnodes = (x, y)
            return {(v, vnodes)}
        else:
            return set()

    def add(self, v, context):
        if len(context) == 0:
            # v is already node, we have to do nothing!
            return
        elif len(context) == 1:
            # replace a node with v
            x, = context
            i = self.nodes.index(x)
            self.nodes[i] = v
            # replace two adjacent edges
            self.edges[i-1] = (self.edges[i-1][0], v)
            self.edges[i] = (v, self.edges[i][1])
        elif len(context) == 2:
            # insert between two nodes
            x, y = context
            i = self.nodes.index(y)
            self.nodes.insert(i, v)
            # replace edge x-y with x-v and insert v-y
            self.edges[i-1] = (x, v)
            self.edges.insert(i, (v, y))
        #self.init_edges()

    def close(self, x, y):
        return sum((t - s)**2 for t, s in zip(x, y)) < self.TOL**2

    def plot(self, *args, **kwargs):
        #plt.figure()
        x0, x1 = zip(*(self.nodes + [self.nodes[0]]))
        plt.plot(x0, x1, *args, **kwargs)

    def clip_from_right(self, a1, b1, c0):
        # return polygon a-b-c-d which is bounded at the left a-b side by self
        # and add intersection nodes to self
        a = self.right_intersection(a1)
        b = self.right_intersection(b1)
        c = (c0, b1)
        d = (c0, a1)
        nodes = self.nrange(a, b, -1) + [c, d]
        pol = Polygon(nodes)
        # useful meta-info about polygon
        pol.a = a
        pol.b = b
        pol.c = c
        pol.d = d
        return pol

    def clip_from_left(self, d1, c1, a0):
        # return polygon a-b-c-d which is bounded at the right c-d side by self
        # and add intersection nodes to self
        c = self.left_intersection(c1)
        d = self.left_intersection(d1)
        a = (a0, d1)
        b = (a0, c1)
        nodes = [a, b] + self.nrange(c, d, -1)
        pol = Polygon(nodes)
        # useful meta-info about polygon
        pol.a = a
        pol.b = b
        pol.c = c
        pol.d = d
        return pol
    
    def clip_from(self, side, a1, b1, c0):
        # generalization of above two
        # side in {"top", "bottom", "left", "right"}
        # with naming logic as in cip_from_right
        # (a1, b1 are clockwise, not ascending)
        a = self.add_intersection(side, a1)
        b = self.add_intersection(side, b1)
        c = (c0, b1) if side in {"left", "right"} else (b1, c0)
        d = (c0, a1) if side in {"left", "right"} else (a1, c0)
        pol = Polygon(self.nrange(a, b, -1) + [c, d])
        a, b, c, d = canonical_order(a, b, c, d)
        pol.set_corners(a, b, c, d)
        return pol

    def cut_path(self, a, b):
        n = self.len()
        ia = self.index(a)
        ib = self.index(b)
        self.edges[ia] = (a, b)
        nodes = list(self.nodes)
        edges = list(self.edges)
        for i in irange(ia, ib, n, 1)[1:-1]:
            #print i
            self.nodes.remove(nodes[i])
            self.edges.remove(edges[i])

    def cut_from_right(self, r):
        # TODO: fails if cut polygon is not connected any more
        nodes = self.all_intersections(r, 0)
        nodes = sorted(nodes, key=lambda x: x[1])[::-1]
        n = self.len()
        for i, node in enumerate(nodes):
            # check whether nodes in between are left or right of line
            j = self.nodes.index(node)
            nextnode = self.nodes[(j+1) % n]
            if nextnode[0] > node[0]:
                self.cut_path(node, nodes[i+1])

    def cut_from_left(self, r):
        nodes = self.all_intersections(r, 0)
        nodes = sorted(nodes, key=lambda x: x[1])
        n = self.len()
        for i, node in enumerate(nodes):
            # check whether nodes in between are left or right of line
            j = self.nodes.index(node)
            nextnode = self.nodes[(j+1) % n]
            if nextnode[0] < node[0]:
                self.cut_path(node, nodes[i+1])

    def split(self, r):
        p, q = self.copy(), self.copy()
        p.cut_from_right(r)
        q.cut_from_left(r)
        return p, q

    def rmin(self):
        return min(self.nodes, key=lambda v: v[0])

    def rmax(self):
        return max(self.nodes, key=lambda v: v[0])

    def zmin(self):
        return min(self.nodes, key=lambda v: v[1])

    def zmax(self):
        return max(self.nodes, key=lambda v: v[1])

    def radiustop(self):
        zmax = max(x[1] for x in self.nodes)
        return min(x[0] for x in self.nodes if x[1]==zmax)

    def radiusbottom(self):
        zmin = min(x[1] for x in self.nodes)
        return min(x[0] for x in self.nodes if x[1]==zmin)

    def set_corners(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def nrange(self, a, b, sign=1):
        # return nodes that range from a to b (including b)
        # sign = 1 for clockwise, -1 for counter-clockwise
        n = self.len()
        ia = self.index(a)
        ib = self.index(b)
        return [self.nodes[i] for i in irange(ia, ib, n, sign)]

    def edgerange(self, a, b):
        # return set of edges from a to b
        nodes = self.nrange(a, b)
        return set(nodes2edges(nodes))

    def index(self, x):
        return self.nodes.index(x)

    def len(self):
        return len(self.nodes)

class MultiPolygon(Polygon):
    "collection of polygons that can be modified together"
    # polygons are list and corresponding outer boundaries are
    # kept in list with same order
    # the Polygon itself is the disjoint union
    def __init__(self, *polygons):
        self.polygons = polygons
        unions, boundaries = compute_disjoint_union(*polygons)
        # at the moment, only support simply connected polygons
        assert len(unions) == 1
        self.boundaries = boundaries
        self.n = len(polygons)
        Polygon.__init__(self, unions[0].nodes)

    def plot(self, *args, **kwargs):
        colors = ["b", "g", "r", "y"]
        colors = [colors[i % 4] for i in range(self.n)]
        for i in range(self.n):
            plot_edges(self.boundaries[i], "-%s" % (colors[i],))
        Polygon.plot(self, *args, **kwargs)

    def add(self, v, context):
        # add node to union and parts; assume union topology does not change
        u = self

        if len(context) == 0:
            # v is already node, we have to do nothing!
            return
        elif len(context) == 1:
            # replace a node with v
            x, = context
            # replace node in polygons containing it
            for p in self.polygons:
                if x in p.nodes:
                    p.add(v, context)
            # replace node in union
            i = u.nodes.index(x)
            e0 = u.edges[i-1]
            e1 = u.edges[i]
            u.nodes[i] = v
            # replace two adjacent boundary edges
            u.edges[i-1] = (u.edges[i-1][0], v)
            u.edges[i] = (v, u.edges[i][1])
            # replace the same edges also in boundaries
            for b in self.boundaries:
                if e0 in b:
                    b.remove(e0)
                    b.add(u.edges[i-1])
                if e1 in b:
                    b.remove(e1)
                    b.add(u.edges[i])
        elif len(context) == 2:
            # insert between two nodes
            # this means new node is on boundary edge, which belongs
            # to exactly one polygon # but let's not use this
            x, y = context
            for p in self.polygons:
                if x in p.nodes and y in p.nodes:
                    p.add(v, context)
            i = u.nodes.index(y)
            u.nodes.insert(i, v)
            # replace edge x-y with x-v and insert v-y
            e0 = u.edges[i-1]
            u.edges[i-1] = (x, v)
            u.edges.insert(i, (v, y))
            for b in self.boundaries:
                if e0 in b:
                    b.remove(e0)
                    b |= {u.edges[i-1], u.edges[i]}

    def all_intersections(self, z, axis=1):
        X = Polygon.all_intersections(self, z, axis)
        for p in self.polygons:
            X1 = p.all_intersections(z, axis)
            X.extend(X1)
        return X

    def cut_from_right(self, r):
        if r <= self.rmax()[0]:
            for p in self.polygons:
                p.cut_from_right(r)
            self.__init__(*self.polygons)

    def cut_from_left(self, r):
        if r <= self.rmax()[0]:
            for p in self.polygons:
                p.cut_from_left(r)
            self.__init__(*self.polygons)

class Ball(object):

    def __init__(self, x0, r, lc=1.):
        self.x0 = x0
        self.r = r
        self.lc = lc
        # TODO: this makes only sense for 2D:
        self.z = z = x0[-1]
        nodes = [(0, z-r), (0, z), (0, z+r)]
        self.nodes = nodes
        x1, x2, x3 = tuple(nodes)
        self.edges = [(x3, x2), (x2, x1), (x1, x2, x3)]

    def boundary(self, dim):
        if dim == 2:
            return {self.edges[2]}
        elif dim == 3:
            return {self}
        elif dim == 1:
            return set()

    def inside(self, x, radius=0.):
        m = np.array([self.x0])
        return np.sum((x - m)**2, 1) <= (self.r + radius)**2

    def inside_single(self, x, radius=0.):
        m = np.array(self.x0)
        return np.sum((x - m)**2) <= (self.r + radius)**2

    def __repr__(self):
        return "Ball(%s, %s, lc=%s)" % (self.x0, self.r, self.lc)

class EmptySet(object):
    def __init__(self):
         self.nodes = []
         self.edges = []

    def boundary(self, *args, **kwargs):
        return set()

    def __repr__(self):
        return "EmptySet()"

def isempty(poly):
    return isinstance(poly, EmptySet)

def nodes2edges(nodes, closed=False):
    if closed:
        return zip(nodes, nodes[1:] + nodes[0:1])
    else:
        return zip(nodes[:-1], nodes[1:])

class PolygonPore(object):
    def __init__(self, poly, name="protein", **params):
        self.name = name
        self.protein = Polygon(poly)
        self.polygons = OrderedDict({name: self.protein})
        self.balls = OrderedDict()
        self.params = Params(params)
        self.boundaries = OrderedDict()

        # cs ... protein crosssections for partition of boundary
        self.add_proteinpartition()
        self.add_molecule_ball()

    def add_molecule_ball(self):
        if "x0" in self.params and self.params.x0 is not None:
            self.molecule = True
            lc = self.params.lcMolecule if "lcMolecule" in self.params else 1.
            molecule = Ball(self.params.x0, self.params.rMolecule, lc)
        else:
            self.molecule = False
            molecule = EmptySet()
        self.balls.update(molecule=molecule)

    def add_balls(self, **balls):
        for pname, p in balls.items():
            if not isinstance(p, Ball):
                x0, r = p
                balls[pname] = Ball(x0, r)
        self.balls.update(balls)

    def build_polygons(self):
        R = self.params.R
        if not "Hbot" in self.params:
            H = self.params.H
            self.params["Htop"] = H/2.
            self.params["Hbot"] = H/2.
        Hbot = self.params.Hbot
        Htop = self.params.Htop

        self.add_membrane()

        sections = self.add_poresections()
        sections = self.add_poreregions(sections)
        self.add_bulkfluids(R, Htop, Hbot, sections)
        for ball in self.balls.values():
            self.add_molecule(ball)

        # for easier handling of alls domains
        self.domains = self.balls.copy()
        self.domains.update(self.polygons)
        return self.polygons

    def build_boundaries(self):
        boundaries = self.boundaries

        # membrane boundary
        mem = self.polygons["membrane"]
        if not isempty(mem):
            memb = mem.edgerange(mem.b, mem.c) | mem.edgerange(mem.d, mem.a)
            boundaries.update(memb=memb)

        # upper, lower, side
        btop = self.polygons["bulkfluid_top"]
        bbot = self.polygons["bulkfluid_bottom"]
        upperb = btop.edgerange(btop.b, btop.c)
        lowerb = bbot.edgerange(bbot.d, bbot.a)
        sideb = btop.edgerange(btop.c, btop.d) | bbot.edgerange(bbot.c, bbot.d)

        boundaries.update(upperb=upperb, lowerb=lowerb, sideb=sideb)
        boundaries.update(self.compute_proteinboundary())

        # molecule
        for bname, ball in self.balls.items():
            boundaries.update({bname + "b": ball.boundary(self.params.dim)})

        return boundaries

    def molecule_intersects(self, z):
        eps = 0.5
        x0 = self.params.x0 if "x0" in self.params else None
        if x0 is None:
            return False
        z0 = x0[-1]
        r = self.params.rMolecule
        return z0 - r - eps <= z <= z0 + r + eps

    def where_is_molecule(self, ball=None):
        "name of domain where ball lies (or None)"
        # TODO: currently this is only based on z position!!
        if ball is None:
            ball = self.balls["molecule"]
        if isempty(ball):
            return None

        domains = ["pore%d" % i for i in range(self.nsections)]
        if self.params.poreregion:
            domains = ["poreregion_bottom"] + domains + ["poreregion_top"]
        else:
            domains = ["bulkfluid_bottom"] + domains + ["bulkfluid_top"]
        i0 = bisect(self.cs, ball.x0[-1])
        return domains[i0]

    def add_molecule(self, ball):
        # do nothing if empty
        if isempty(ball):
            return

        # check in which domain ball is contained
        domstr = self.where_is_molecule(ball)
        domain = self.polygons[domstr]

        # in 3D, add ball as hole in containing domain
        if "dim" in self.params and self.params.dim == 3:
            if not hasattr(domain, "holes"):
                domain.holes = []
            domain.holes.append(ball)
            return

        # in 2D, modify containing domain
        x1, x2, x3 = tuple(ball.nodes)
        domain.left_intersection(x1[1])
        domain.left_intersection(x3[1])
        i = domain.edges.index((x1, x3))
        domain.edges[i] = (x1, x2, x3)

    def add_membrane(self):
        R = self.params.R
        if not ("no_membrane" in self.params and self.params.no_membrane) and (
                self.protein.rmax()[0] < R):
            hmem = self.params.hmem
            zmem = self.params.zmem
            zbot = zmem - 0.5*hmem
            ztop = zmem + 0.5*hmem
            self.membrane = self.protein.clip_from_right(zbot, ztop, R)
        else:
            self.membrane = EmptySet()
        self.polygons["membrane"] = self.membrane
        return self.membrane

    def add_poresections(self, remove_intersections=True):
        # cs ... crosssections
        cs = self.params.cs if "cs" in self.params else None
        ztop = max(x[1] for x in self.protein.nodes)
        zbot = min(x[1] for x in self.protein.nodes)
        if cs is None:
            cs = [zbot + i/3.*(ztop - zbot) for i in [1, 2]]
        else:
            cs = list(cs)
        assert all(zbot < z < ztop for z in cs)
        cs = [zbot] + sorted(cs) + [ztop]
        # do not add crosssections that intersect with molecule
        if remove_intersections:
            for z in list(cs):
                if self.molecule_intersects(z):
                    cs.remove(z)

        pairs = zip(cs[:-1], cs[1:])
        sections = tuple([self.protein.clip_from_left(a, b, 0) for a, b in pairs])
        names = ["pore%d" % i for i in range(len(sections))]
        self.nsections = len(sections)
        self.lpore = ztop - zbot
        self.cs = cs
        self.polygons.update({name: s for name, s in zip(names, sections)})
        return sections

    def add_poreregions(self, sections):
        if not "poreregion" in self.params or not self.params.poreregion:
            self.params.poreregion = False
            self.polygons["poreregion_top"] = EmptySet()
            self.polygons["poreregion_bottom"] = EmptySet()
            return sections
        if not "H0" in self.params or self.params.H0 is None:
            self.params["H0"] = 0.5*self.params.H + self.protein.zmax()[1]
        if not "R0" in self.params or self.params.R0 is None:
            self.params["R0"] = self.protein.rmax()[0]
        H0 = self.params.H0
        R0 = self.params.R0
        H0top = H0bot = H0/2.
        #print self.protein.zmax()[1], H0top, self.params.H*0.5

        center = MultiPolygon(self.protein, self.membrane, *sections)
        prtop = center.clip_from("top", R0, 0, H0top)
        prbot = center.clip_from("bottom", 0, R0, -H0bot)

        # OLD VERSION: doesn't work bc pore region and memb can overlap
#        section = sections[-1]
#        a = section.b
#        d = self.protein.top_intersection(R0)
#        upper = self.protein.nrange(d, section.c, -1) + [a]
#        b, c = (0, H0top), (R0, H0top)
#        nodes = [b, c] + upper
#        prtop = Polygon(nodes)
#        prtop.set_corners(a, b, c, d)
#
#        section = sections[0]
#        b = section.a
#        c = self.protein.bottom_intersection(R0)
#        lower = [b] + self.protein.nrange(section.d, c, -1)
#        a, d = (0, -H0bot), (R0, -H0bot)
#        nodes = lower + [d, a]
#        prbot = Polygon(nodes)
#        prbot.set_corners(a, b, c, d)

        self.polygons["poreregion_top"] = prtop
        self.polygons["poreregion_bottom"] = prbot
        sections = [prbot] + list(sections) + [prtop]
        return sections

    def add_bulkfluids(self, R, Htop, Hbot, sections):
        # sections ... pore sections that remain after adding molecule
        polygons = list(sections) + [self.protein, self.membrane]

        upper = compute_upper_boundary(*polygons)
        b, c = (0, Htop), (R, Htop)
        nodes = [b, c] + upper
        btop = Polygon(nodes)
        btop.a = upper[-1]
        btop.b = b
        btop.c = c
        btop.d = upper[0]

        lower = compute_lower_boundary(*polygons)
        d, a = (R, -Hbot), (0, -Hbot)
        nodes = lower + [d, a]
        bbot = Polygon(nodes)
        bbot.a = a
        bbot.b = lower[0]
        bbot.c = lower[-1]
        bbot.d = d

        self.polygons["bulkfluid_top"] = btop
        self.polygons["bulkfluid_bottom"] = bbot
        return btop, bbot

    def add_proteinpartition(self):
        # do at the beginning
        # insert intersection points with cs in polygon
        cs = self.params.proteincs if "proteincs" in self.params else None
        if cs is None or len(cs)==0:
            self.proteinpartition = False
            return
        protein = self.protein
        ztop = max(x[1] for x in protein.nodes)
        zbot = min(x[1] for x in protein.nodes)
        assert all(zbot < z < ztop for z in cs)
        for z in cs:
            protein.all_intersections(z)
        self.proteincs = cs
        self.proteinpartition = True

    def compute_proteinboundary(self):
        # do after modifying edges!!
        protein = self.protein
        if not self.proteinpartition:
            dic = {"%sb" % self.name: set(self.protein.edges)}

        else:
            cs = self.proteincs
            ztop = max(x[1] for x in protein.nodes)
            zbot = min(x[1] for x in protein.nodes)
            cs = [zbot] + sorted(cs) + [ztop]
            npart = len(cs) - 1
            sets = [set() for i in range(npart)]

            # partition edges
            for edge in protein.edges:
                x = min(edge, key=lambda t: t[1])
                ix = bisect(cs, x[1]) - 1
                if ix == npart:
                    ix -= 1
                sets[ix].add(edge)

            dic = {"%sb%d" % (self.name, i): sets[i] for i in range(npart)}

        # subtract edges on protein-membrane interface
        mem = self.polygons["membrane"]
        if not isempty(mem):
            pmemb = protein.edgerange(mem.b, mem.a)
            for s in dic.values():
                s -= pmemb

        return dic
    
    def radius_at(self, z):
        "return z-dependent pore radius"
        I = self.protein.intersections(z)
        return min(x[0] for x in I.keys())

class MultiPolygonPore(PolygonPore):
    """For pores consisting of multiple subdomains, e.g. Wei-Rant pore,
    but without further pore bisection (proteinpartition).
    Assumes that polygons are disjoint and their union simply connected."""
    def __init__(self, **params):
        self.polygons = OrderedDict()
        self.balls = OrderedDict()
        self.params = Params(params)
        self.boundaries = OrderedDict()
        self.names = []
        self.add_molecule_ball()

    def add_polygons(self, **polygons):
        "enables adding polygons in defined order"
        for pname, p in polygons.items():
            if not isinstance(p, Polygon):
                polygons[pname] = Polygon(p)
        self.polygons.update(polygons)
        self.names.extend(polygons.keys())

    def build_polygons(self):
        # read global height, width params
        R = self.params.R
        if not "Hbot" in self.params:
            H = self.params.H
            self.params["Htop"] = H/2.
            self.params["Hbot"] = H/2.
        Hbot = self.params.Hbot
        Htop = self.params.Htop

        # first build encompassing polygon out of existing.
        self.protein = MultiPolygon(*self.polygons.values())
        self.protein.cut_from_right(R)
        self.polygons = OrderedDict()

        bnames = ["%sb" % name for name in self.names]
        self.proteinb = OrderedDict(zip(bnames, self.protein.boundaries))

        self.add_membrane()
        sections = self.add_poresections()
        sections = self.add_poreregions(sections)
        self.add_bulkfluids(R, Htop, Hbot, sections)
        self.polygons.update(zip(self.names, self.protein.polygons))
        for ball in self.balls.values():
            self.add_molecule(ball)

        # for easier handling of alls domains
        self.domains = self.balls.copy()
        self.domains.update(self.polygons)
        return self.polygons

    def build_boundaries(self):
        boundaries = self.boundaries
        proteinb = self.proteinb

        # upper, lower, side
        btop = self.polygons["bulkfluid_top"]
        bbot = self.polygons["bulkfluid_bottom"]
        upperb = btop.edgerange(btop.b, btop.c)
        lowerb = bbot.edgerange(bbot.d, bbot.a)
        sideb = btop.edgerange(btop.c, btop.d) | bbot.edgerange(bbot.c, bbot.d)
        boundaries.update(upperb=upperb, lowerb=lowerb, sideb=sideb)

        fluid_polys = [btop, bbot] + \
            [self.polygons["poreregion_top"],
             self.polygons["poreregion_bottom"]] + \
            [self.polygons["pore%d" % i] for i in range(self.nsections)]
        fluid_edges = set([e[::-1] for p in fluid_polys for e in p.edges])

        # membrane-fluid boundary
        mem = self.polygons["membrane"]
        if not isempty(mem):
            memb = fluid_edges & (mem.edgerange(mem.b, mem.c) | mem.edgerange(mem.d, mem.a))
            boundaries.update(memb=memb)
        else:
            boundaries.update(memb=set())
        # protein-fluid boundary, divided between protein parts
        for bname in proteinb:
            proteinb[bname] &= fluid_edges
        boundaries.update(proteinb)

        # molecule
        for bname, ball in self.balls.items():
            boundaries.update({bname + "b": ball.boundary(self.params.dim)})

        return boundaries

def canonical_order(*corners):
    """return corners of enclosing rectangle in their canonical order 
    (clockwise starting from bottom left).
    can be used to simply reorder rectangle nodes."""
    xmin = min(x[0] for x in corners)
    xmax = max(x[0] for x in corners)
    ymin = min(x[1] for x in corners)
    ymax = max(x[1] for x in corners)
    return (xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)

def join_nodes(polygons):
    return list(set([x for p in polygons for x in p.nodes]))

def join_edges(polygons):
    return list(set([x for p in polygons for x in p.edges]))

def joint_boundary(polygons):
    "determine edges that do not feature twice"
    edges = set()
    for p in polygons:
        for e in p.edges:
            if e[::-1] in edges:
                edges.remove(e[::-1])
            else:
                edges.add(e)
    return list(edges)

def top_edge_from_node(node, edges):
    edges = [e for e in edges if e[0]==node]
    assert len(edges) > 0
    return max(edges, key=lambda e: e[1][1])

def bottom_edge_from_node(node, edges):
    edges = [e for e in edges if e[1]==node]
    assert len(edges) > 0
    return min(edges, key=lambda e: e[0][1])

def compute_upper_boundary(*polygons):
    allnodes = join_nodes(polygons)
    edges = joint_boundary(polygons)
    # start at top node with r=0, go until rmax
    x0 = max([x for x in allnodes if x[0]==0], key=lambda x: x[1])
    rmax = max([x[0] for x in allnodes])
    X = [x0]
    while x0[0] < rmax:
        edge = top_edge_from_node(x0, edges)
        x0 = edge[1]
        X.append(x0)
    return X[::-1]

def compute_lower_boundary(*polygons):
    allnodes = join_nodes(polygons)
    edges = joint_boundary(polygons)
    # start at bottom node with r=0, go until rmax
    x0 = min([x for x in allnodes if x[0]==0], key=lambda x: x[1])
    rmax = max([x[0] for x in allnodes])
    X = [x0]
    while x0[0] < rmax:
        edge = bottom_edge_from_node(x0, edges)
        x0 = edge[0]
        X.append(x0)
    return X

def compute_disjoint_union(*polygons):
    """compute boundary of disjoint union of polygons,
    return union polygon(s) and list of boundary pieces in the same order as
    corresponding input polygons"""
    # note: the outputs that are holes are correctly returned in clockwise
    # i.e. opposite direction. this can be identified with the is_a_hole
    # function below.

    # determine edges that do not feature twice
    edges = {} # boundary edge: polygon index
    boundaries = [set() for i in polygons]
    for i, ply in enumerate(polygons):
        for e in ply.edges:
            e1 = e[::-1]
            if e1 in edges:
                i1 = edges[e1]
                edges.pop(e1)
                boundaries[i1].remove(e1)
            else:
                edges[e] = i
                boundaries[i].add(e)

    polys = []
    # construct connected subsets of unit polygon
    # start with any boundary edge and iteratively add all connected edges
    # repeat until no edges left
    while edges:
        e = next(iter(edges))
        edges.pop(e)
        v0 = e[0]
        v = e[-1]
        nodes = [v0]
        while v != v0:
            nodes.append(v)
            e = next(e1 for e1 in edges if e1[0]==v)
            edges.pop(e)
            v = e[-1]
        polys.append(Polygon(nodes))

    return polys, boundaries

def is_a_hole(poly):
    """determine whether polygon is a hole (edges run counter-clockwise) or not
    by using the showlace formula"""
    # s = twice the signed area
    s = sum((e[-1][0] - e[0][0])*(e[-1][1] + e[0][1]) for e in poly.edges)
    return s < 0.

"""winding number algorithm from
http://geomalgorithms.com/a03-_inclusion.html,
implemented for array of points at once."""

def winding_number(poly, x):
    "return array of winding numbers of poly around points x"
    wn = np.zeros(x.shape[0], dtype=int)
    for e in poly.edges:
        v0, v1 = e
        lor = left_on_right(v0, v1, x)
        is_left = lor > 0
        is_right = lor < 0
        upward_crossing = (v0[1] <= x[:, 1]) & (v1[1] > x[:, 1])
        downward_crossing = (v0[1] > x[:, 1]) & (v1[1] <= x[:, 1])
        wn[upward_crossing & is_left] += 1
        wn[downward_crossing & is_right] -= 1
    return -wn # minus because winding convention of edges is wrong

def left_on_right(v0, v1, x):
    """tests whether point x is left (output>0), on (=0), or right (<0)
    of the infinite line defined by v0, v1. x is many points as (N, 2) array"""
    return (v1[0] - v0[0])*(x[:,1] - v0[1]) - (x[:,0] - v0[0])*(v1[1] - v0[1])

#def is_left(v0, v1, x):
#    return left_on_right(e, x) > 0.
#
#def is_right(v0, v1, x):
#    return left_on_right(e, x) < 0.

def union_example():
    A = Polygon([(0., 0.), (0., 3.), (1., 3.), (1., 2.), (1., 1.), (1., 0.)])
    B = Polygon([(1., 2.), (1., 3.), (2., 3.), (2., 2.)])
    C = Polygon([(1., 0.), (1., 1.), (2., 1.), (2., 0.)])
    D = Polygon([(2., 0.), (2., 1.), (2., 2.), (2., 3.), (3., 3.), (3., 0.)])
#    unions, b = compute_disjoint_union(A, B, C, D)
#
#    for i, c in enumerate(["g", "b", "r", "y"]):
#        plot_edges(b[i], "-%s" % c)
#    unions[0].plot(".k")
#    unions[1].plot("og")
#    from matplotlib.pyplot import xlim, ylim
#    xlim(-1, 4)
#    ylim(-1, 4)
    return A, B, C, D


def plot_edges(edges, *args, **kwargs):
    for t in edges:
        if len(t) == 3:
            x, _, y = t
        else:
            x, y = t
        plt.plot([x[0], y[0]], [x[1], y[1]], *args, **kwargs)

if __name__ == "__main__":
    from nanopores.geometries.alphahempoly import poly
    params = dict(
        R = 4.1,
        H = 30,
        hmem = 2.2,
        zmem = -6,
        cs = [-3.3, -6.6],
        proteincs=[-2.5, -4.9, -7.51],
        x0 = [0.,0.,-6],
        rMolecule = .5,
        dim = 2,
        #no_membrane = True
    )

    #p = PolygonPore(poly, "ahem", **params)
    p = MultiPolygonPore(**params)
    p.add_polygons(ahem=poly)
    p.build_polygons()
    p.build_boundaries()
    print p.polygons.keys()
    print p.boundaries.keys()

    print "mol in", p.where_is_molecule()
    p.protein.plot(".k", zorder=100)
    p.polygons["bulkfluid_top"].plot("--r")
    p.polygons["bulkfluid_bottom"].plot("--g")
    p.polygons["pore0"].plot(".-b")

    plot_edges(p.boundaries["memb"], color="b")
    plot_edges(p.boundaries["lowerb"], color="r")
    plot_edges(p.boundaries["upperb"], color="g")
    plot_edges(p.boundaries["sideb"], color="y")
    plt.xlim(-0.2, params["R"] + 0.2)
    plt.ylim(-11, 1)

    #plot_edges(p.boundaries["ahemb0"], color="red")
    #plot_edges(p.boundaries["ahemb1"], color="yellow")
    #plot_edges(p.boundaries["ahemb2"], color="green")
    #plot_edges(p.boundaries["ahemb3"], color="red")
    plt.show()
#    r = 6.
#    P = p.protein
#    P.plot(".--k")
#    plt.axvline(x=r)
#    plt.show()
#    P.all_intersections(r, 0)
#    P.plot(".--k")
#    plt.axvline(x=r)
#    plt.show()
#    #print [v for v in P.nodes if v[0]==r]
#    P.cut_from_right(r)
#    P.plot(".--k")
    #P.cut_from_right(r)
    #P.plot(".-")
