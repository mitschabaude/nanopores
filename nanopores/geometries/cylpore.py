# (c) 2017 Gregor Mitscha-Baude
import nanopores.py4gmsh as gmsh
from nanopores.geo2xml import geofile2geo
from nanopores.tools.polygons import PolygonPore, plot_edges, isempty
from nanopores.tools.utilities import Log

default = dict(
    geoname = "cylpore",
    porename = "protein",
    dim = 2,
    R = 15.,
    H = 30.,
    x0 = None,
    rMolecule = 0.5,
    lcMolecule = 0.25,
    lcCenter = 0.5,
    hmem = 2.2,
    zmem = 0.,
    cs = (), # crosssections; list of z coordinates
    poreregion = False, # whether to include fluid above pore as subdomain
)

class Pore(PolygonPore):
    "PolygonPore with additional gmsh creation functions"
    def __init__(self, poly, **params):
        params = dict(default, **params)
        self.dim = params["dim"]
        self.geoname = params["geoname"] # name of folder where mesh is created
        name = params["porename"] # name of pore wall subdomain, i.e. "dna"
        PolygonPore.__init__(self, poly, name, **params)

    def build(self, h=1.):
        with Log("computing polygons..."):
            self.build_polygons()
            self.build_boundaries()
        geo = self.build_geometry(h)
        return geo

    def build_geometry(self, h=1.):
        with Log("writing gmsh code..."):
            code, meta = self.to_gmsh(h)
        geo = geofile2geo(code, meta, name=self.geoname)
        return geo

    def to_gmsh(self, h=1.):
        # create lists with common nodes, edges
        # TODO: maybe we should move away from indices and implement map
        # {entity: gmsh_entity} ???
        # set lookup is probably faster than list lookup !!!
        self.nodes = list(set([x for p in self.polygons.values() for x in p.nodes]))
        self.edges = list(set([e for p in self.polygons.values() for e in p.edges]))
        self.gmsh_nodes = [None for x in self.nodes]
        self.gmsh_edges = [None for x in self.edges]

        # in 2D, create plane surfaces for every polygon
        # TODO: 3D
        dim = self.dim
        lcs = self.set_length_scales(h)

        for pname, p in self.polygons.items():
            gmsh.Comment("Creating %s subdomain." %pname)
            if isempty(p):
                gmsh.NoPhysicalVolume(pname)
                continue
            if dim == 2:
                lc = lcs[pname]
                #print pname, lc
                ll = self.LineLoop(p.edges, lc)
                vol = gmsh.PlaneSurface(ll)
                gmsh.PhysicalVolume(vol, pname, dim)
            else:
                raise NotImplementedError

        # add physical surfaces
        for bname, bset in self.boundaries.items():
            self.PhysicalBoundary(bset, bname)

        gmsh.raw_code(["General.ExpertMode = 1;"])
        #gmsh.raw_code(["Mesh.Algorithm3D = 2;"])
        code = gmsh.get_code()
        meta = gmsh.get_meta()
        self.gmshcode = code, meta
        return code, meta

    def set_length_scales(self, h):
        lc = {p: h for p in self.polygons}
        lc["molecule"] = h*self.params.lcMolecule
        for i in range(self.nsections):
            lc["pore%d" % i] = h*self.params.lcCenter
        else:
            lc["%s" % (self.name,)] = h*self.params.lcCenter
        return lc

    def Point(self, x, lc):
        i = self.nodes.index(x)
        gmsh_x = self.gmsh_nodes[i]
        if gmsh_x is not None:
            return gmsh_x

        #print x, lc
        x = x + tuple(0. for i in range(3 - self.dim))
        gmsh_x = gmsh.Point(x, lc)
        self.gmsh_nodes[i] = gmsh_x
        return gmsh_x

    def Edge(self, e, lc=1.):
        # generalizes Line and Circle
        # if exists, return gmsh edge
        i = self.edges.index(e)
        gmsh_e = self.gmsh_edges[i]
        if gmsh_e is not None:
            return gmsh_e
        # otherwise, discriminate between Line and Circle
        points = [self.Point(v, lc) for v in e]

        if len(e) == 2:
            gmsh_e = gmsh.Line(*points)
        elif len(e) == 3:
            gmsh_e = gmsh.Circle(points)

        self.gmsh_edges[i] = gmsh_e
        # also save flipped edge (if exists)
        try:
            i1 = self.edges.index(e[::-1])
            self.gmsh_edges[i1] = "-" + gmsh_e
        except ValueError:
            pass
        return gmsh_e

    def LineLoop(self, ll, lc=1.):
        lines = [self.Edge(edge, lc) for edge in ll]
        ll = gmsh.LineLoop(lines)
        return ll

    def PhysicalBoundary(self, bset, bname):
        gmsh.Comment("Creating %s boundary." %bname)
        dim = self.dim
        if not bset:
            gmsh.NoPhysicalSurface(bname)
        if dim == 2:
            boundary = [self.Edge(e) for e in bset]
            gmsh.PhysicalSurface(boundary, bname, dim)
        else:
            raise NotImplementedError


def get_geo(poly, h=1., **params):
    p = Pore(poly, **params)
    geo = p.build(h=h)
    return geo

if __name__ == "__main__":
    from alphahempoly import poly
    from matplotlib import pyplot as plt
    from dolfin import plot, interactive
    from nanopores import user_params
    params = user_params(
        h = 1.,
        porename="ahem",
        cs=[-3, -6],
        zmem=-5.,
        proteincs=[-5.],
        #x0 = [0., 0., 0.],
    )

    geo = get_geo(poly, **params)
    plot(geo.subdomains)
    plot(geo.boundaries)
    print geo
    interactive()

    # --- TEST
    p = Pore(poly, **params)
    p.build_polygons()
    p.build_boundaries()
    print p.polygons.keys()
    print p.boundaries.keys()

    p.protein.plot(".k", zorder=100)
    #p.polygons["bulkfluid_top"].plot()
    p.polygons["pore0"].plot()

    plot_edges(p.boundaries["memb"], color="blue")
    plot_edges(p.boundaries["lowerb"])
    plot_edges(p.boundaries["upperb"])
    plot_edges(p.boundaries["sideb"], color="yellow")
    plot_edges(p.boundaries["ahemb0"], color="red")
    #plt.show()
    # ---