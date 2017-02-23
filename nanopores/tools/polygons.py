# (c) 2017 Gregor Mitscha-Baude
"""this was intended to perform boolean operations on polygons,
but this seems to be too hard without relying on external software"""
from matplotlib import pyplot as plt
#from nanopores.tools.utilities import union, collect

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
            v = (x0 + (y0 - x0)*t, x1 + (y1 - x1)*t)
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
            # v is already node, we have to do nothing!"o-k"
            return
        elif len(context) == 1:
            # replace a node with v
            x, = context
            i = self.nodes.index(x)
            self.nodes[i] = v
        elif len(context) == 2:
            # insert between two nodes
            x, y = context
            i = self.nodes.index(y)
            self.nodes.insert(i, v)
        self.init_edges()

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
        return Polygon(nodes)

    def clip_from_left(self, d1, c1, a0):
        # return polygon a-b-c-d which is bounded at the right c-d side by self
        # and add intersection nodes to self
        c = self.left_intersection(c1)
        d = self.left_intersection(d1)
        a = (a0, d1)
        b = (a0, c1)
        nodes = [a, b] + self.nrange(c, d, -1)
        return Polygon(nodes)

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

    def index(self, x):
        return self.nodes.index(x)

    def len(self):
        return len(self.nodes)

class PorePolygon(Polygon):
    def get_membrane(self, hmem, zmem, R):
        zbot = zmem - 0.5*hmem
        ztop = zmem + 0.5*hmem
        return self.clip_from_right(zbot, ztop, R)

    def get_poresections(self, cs=None):
        # cs ... crosssections
        ztop = max(x[1] for x in self.nodes)
        zbot = min(x[1] for x in self.nodes)
        cs = [] if cs is None else list(cs)
        assert all(zbot < z < ztop for z in cs)
        cs = [zbot] + sorted(cs) + [ztop]
        pairs = zip(cs[:-1], cs[1:])
        return tuple([self.clip_from_left(a, b, 0) for a, b in pairs])

    def get_bulkfluids(self, R, H, *polygons):
        upper = compute_upper_boundary(*polygons)
        lower = compute_lower_boundary(*polygons)
        nodes = [(0, H/2.), (R, H/2.)] + upper
        bulkfluid_top = Polygon(nodes)
        nodes = lower + [(R, -H/2.), (0, -H/2.)]
        bulkfluid_bottom = Polygon(nodes)
        return bulkfluid_top, bulkfluid_bottom


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

if __name__ == "__main__":
    from nanopores.geometries.alphahempoly import poly
    R = 10
    H = 30
    hmem = 2
    zmem = -6
    cs = [-2, -4]

    p = PorePolygon(poly)
    pmem = p.get_membrane(hmem, zmem, R)
    pbot, pctr, ptop = p.get_poresections(cs)

    p.plot()
    pmem.plot()
    top, bottom = p.get_bulkfluids(R, H, p, pmem, pctr, pbot)
    top.plot()
    bottom.plot()
    plt.show()

