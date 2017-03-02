# (c) 2017 Gregor Mitscha-Baude
"""this was intended to perform boolean operations on polygons,
but this seems to be too hard without relying on external software.
now it is simply for preparing the polygons involved in a pore geometry where
only the pore protein crosssection and parameters have to be plugged in."""
from bisect import bisect
from collections import OrderedDict
from matplotlib import pyplot as plt
from nanopores.tools.utilities import Params

class Polygon(object):
    TOL = 0.1

    def __init__(self, nodes):
        self.nodes = [(float(x[0]), float(x[1])) for x in nodes]
        self.init_edges()

    def init_edges(self):
        nodes = self.nodes
        self.edges = zip(nodes, nodes[1:] + nodes[0:1])

    def intersections(self, z):
        z = float(z)
        return {v: vnodes for edge in self.edges \
                          for v, vnodes in self.intersect_edge(edge, z)}

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

    def all_intersections(self, z):
        I = self.intersections(z)
        for x in I:
            self.add(x, I[x])
        return I.keys()

    def intersect_edge(self, edge, z):
        # return intersection points with line z=z
        # and vertices to provide context for insertion
        x, y = edge
        x0, x1 = x
        y0, y1 = y
        # special case
        if x1 == y1 == z:
            return {(x, ()), (y, ())}
        elif x1 <= z < y1 or y1 < z <= x1:
            # special case
            if abs(y1 - x1) < 1e-3*self.TOL:
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

    def nrange(self, a, b, sign=1):
        # return nodes that range from a to b (including b)
        # sign = 1 for counter-clockwise, -1 for clockwise
        n = self.len()
        ia = self.index(a)
        ib = self.index(b)
        if sign == -1:
            ib = (ib - 1) % n
            if ib > ia:
                ia += n
            ran = range(ia, ib, -1)
        elif sign == 1:
            ib = (ib + 1) % n
            if ia > ib:
                ib += n
            ran = range(ia, ib)
        else:
            raise NotImplementedError
        return [self.nodes[i % n] for i in ran]

    def edgerange(self, a, b):
        # return set of edges from a to b
        nodes = self.nrange(a, b)
        return set(nodes2edges(nodes))

    def index(self, x):
        return self.nodes.index(x)

    def len(self):
        return len(self.nodes)

class HalfCircle(object):
    # could also be implemented with two circular arcs
    def __init__(self, z, r):
        nodes = [(0, z-r), (0, z), (0, z+r)]
        self.nodes = nodes
        x1, x2, x3 = tuple(nodes)
        self.edges = [(x3, x2), (x2, x1), (x1, x2, x3)]

    def boundary(self):
        return {self.edges[2]}

    # TODO: def plot(self):

class EmptySet(object):
    def __init__(self):
         self.nodes = []
         self.edges = []
    def boundary(self):
        return set()

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
        self.params = Params(params)
        self.boundaries = OrderedDict()

        # cs ... protein crosssections for partition of boundary
        self.add_proteinpartition()

        if "x0" in self.params and self.params.x0 is not None:
            self.molecule = True
        else:
            self.molecule = False

    def build_polygons(self):
        hmem = self.params.hmem
        zmem = self.params.zmem
        R = self.params.R
        if "Hbot" in self.params:
            Hbot = self.params.Hbot
            Htop = self.params.Htop
        else:
            H = self.params.H
            Htop = Hbot = H/2.

        self.add_membrane(hmem, zmem, R)

        cs = self.params.cs if "cs" in self.params else []
        sections = self.add_poresections(cs=cs)
        self.add_bulkfluids(R, Htop, Hbot, sections)
        self.add_molecule()

        # TODO: maybe this could be done prettier
        # rebuild polygon dict with new order (molecule first)
        # order matters because some first length scale overrides later ones
        poly = self.polygons
        mol = poly.pop("molecule")
        self.polygons = OrderedDict(molecule=mol)
        self.polygons.update(poly)

        return self.polygons

    def build_boundaries(self):
        boundaries = self.boundaries

        # membrane boundary
        mem = self.polygons["membrane"]
        memb = mem.edgerange(mem.b, mem.c) | mem.edgerange(mem.d, mem.a)

        # upper, lower, side
        btop = self.polygons["bulkfluid_top"]
        bbot = self.polygons["bulkfluid_bottom"]
        upperb = btop.edgerange(btop.b, btop.c)
        lowerb = bbot.edgerange(bbot.d, bbot.a)
        sideb = btop.edgerange(btop.c, btop.d) | bbot.edgerange(bbot.c, bbot.d)

        boundaries.update(memb=memb, upperb=upperb, lowerb=lowerb, sideb=sideb)
        boundaries.update(self.compute_proteinboundary())

        # molecule
        moleculeb = self.polygons["molecule"].boundary()
        boundaries.update(moleculeb=moleculeb)

        return boundaries

    def molecule_intersects(self, z):
        eps = 0.5
        x0 = self.params.x0 if "x0" in self.params else None
        if x0 is None:
            return False
        z0 = x0[-1]
        r = self.params.rMolecule
        return z0 - r - eps <= z <= z0 + r + eps

    def where_is_molecule(self):
        # TODO: currently this is only based on z position!!
        if not self.molecule:
            return None

        domains = ["pore%d" % i for i in range(self.nsections)]
        domains = ["bulkfluid_bottom"] + domains + ["bulkfluid_top"]
        i0 = bisect(self.cs, self.params.x0[-1])
        return domains[i0]

    def add_molecule(self):
        if not self.molecule:
            self.polygons["molecule"] = EmptySet()
            return

        # add molecule to polygons
        z = self.params.x0[-1]
        r = self.params.rMolecule
        molecule = HalfCircle(z, r)
        self.polygons["molecule"] = molecule

        # modify domain in which molecule is contained
        domstr = self.where_is_molecule()
        domain = self.polygons[domstr]
        a = domain.a
        b = domain.b
        x1, x2, x3 = tuple(molecule.nodes)
        domain.add(x1, (a, b))
        domain.add(x3, (x1, b))
        i = domain.edges.index((x1, x3))
        domain.edges[i] = (x1, x2, x3)

    def add_membrane(self, hmem, zmem, R):
        zbot = zmem - 0.5*hmem
        ztop = zmem + 0.5*hmem
        self.membrane = self.protein.clip_from_right(zbot, ztop, R)
        self.polygons["membrane"] = self.membrane
        return self.membrane

    def add_poresections(self, cs=None, remove_intersections=True):
        # cs ... crosssections
        ztop = max(x[1] for x in self.protein.nodes)
        zbot = min(x[1] for x in self.protein.nodes)
        cs = [] if cs is None else list(cs)
        assert all(zbot < z < ztop for z in cs)
        cs = [zbot] + sorted(cs) + [ztop]
        # do not add crosssections that intersect with molecule
        if remove_intersections:
            for z in list(cs):
                if self.molecule_intersects(z):
                    cs.remove(z)

        pairs = zip(cs[:-1], cs[1:])
        sections = tuple([self.protein.clip_from_left(a, b, 0) for a, b in pairs])
#        if len(sections) == 3:
#            names = ["porebot", "porectr", "poretop"]
#        elif len(sections) == 2:
#            names = ["porebot", "poretop"]
#        elif len(sections) == 1:
#            names = ["pore"]
#        else:
        names = ["pore%d" % i for i in range(len(sections))]
        self.nsections = len(sections)
        self.params["lpore"] = ztop - zbot
        self.cs = cs
        self.polygons.update({name: s for name, s in zip(names, sections)})
        return sections

    def add_bulkfluids(self, R, Htop, Hbot, polygons):
        # polygons ... pore sections that remain after adding moleculetop
        polygons = list(polygons) + [self.protein, self.membrane]

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
        pmemb = protein.edgerange(mem.b, mem.a)
        for s in dic.values():
            s -= pmemb

        return dic

def join_nodes(polygons):
    return list(set([x for p in polygons for x in p.nodes]))

def join_edges(polygons):
    return list(set([x for p in polygons for x in p.edges]))

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
    edges = join_edges(polygons)
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
    edges = join_edges(polygons)
    # start at bottom node with r=0, go until rmax
    x0 = min([x for x in allnodes if x[0]==0], key=lambda x: x[1])
    rmax = max([x[0] for x in allnodes])
    X = [x0]
    while x0[0] < rmax:
        edge = bottom_edge_from_node(x0, edges)
        x0 = edge[0]
        X.append(x0)
    return X

def plot_edges(edges, *args, **kwargs):
    for x, y in edges:
        plt.plot([x[0], y[0]], [x[1], y[1]], *args, **kwargs)

if __name__ == "__main__":
    from nanopores.geometries.alphahempoly import poly
    params = dict(
        R = 10,
        H = 30,
        hmem = 2,
        zmem = -6,
        cs = [-3.3, -6.6],
        proteincs=[-6.8],
        x0 = [0.,0.,-6],
        rMolecule = .5,
    )

    p = PolygonPore(poly, "ahem", **params)
    p.build_polygons()
    p.build_boundaries()
    print p.polygons.keys()
    print p.boundaries.keys()

    print "mol in", p.where_is_molecule()
    p.protein.plot(".k", zorder=100)
    #p.polygons["bulkfluid_top"].plot()
    p.polygons["pore0"].plot(".-b")

    plot_edges(p.boundaries["memb"], color="blue")
    plot_edges(p.boundaries["lowerb"])
    plot_edges(p.boundaries["upperb"])
    plot_edges(p.boundaries["sideb"], color="yellow")
    plot_edges(p.boundaries["ahemb"], color="red")
    #plt.show()

