# (c) 2017 Gregor Mitscha-Baude
import nanopores.py4gmsh as gmsh
import nanopores.geometries.curved as curved
from nanopores.geo2xml import geofile2geo, reconstructgeo
from nanopores.tools.polygons import (Ball, Polygon,
                                      PolygonPore, MultiPolygonPore,
                                      plot_edges, isempty)
from nanopores.tools.utilities import Log, union

default = dict(
    geoname = "cylpore",
    porename = "protein",
    dim = 2,
    R = 15.,
    H = 15.,
    x0 = None,
    rMolecule = 0.5,
    lcMolecule = 0.25,
    lcCenter = 0.5,
    hmem = 2.2,
    zmem = 0.,
    # TODO: better default behaviour of crosssections
    cs = (), # crosssections; list of z coordinates
    # maybe TODO: add poreregion (and solve stokes only there?)
    poreregion = False, # whether to include fluid above pore as subdomain

)

default_synonymes = dict(
    #subdomains
    bulkfluid = {"bulkfluid_top", "bulkfluid_bottom"},
    fluid = {"bulkfluid", "pore"},
    solid = {"membrane", "poresolid", "molecule"},
    ions = "fluid",
    #boundaries
    noslip = {"poresolidb", "memb", "moleculeb"},
    bV = "lowerb",
    ground = "upperb",
    bulk = {"lowerb", "upperb"},
    nopressure = "upperb",
)

class Pore(PolygonPore):
    "PolygonPore with additional gmsh creation functions"

    def __init__(self, poly, **params):
        params = dict(default, **params)
        self.dim = params["dim"]
        self.geoname = params["geoname"] # name of folder where mesh is created
        name = params["porename"] # name of pore wall subdomain, i.e. "dna"
        self.names = [name]
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
        meta["params"] = dict(self.params)
        self.geo = geofile2geo(code, meta, name=self.geoname)
        self.finalize_geo(self.geo)
        return self.geo

    def finalize_geo(self, geo):
        geo.params = self.params
        geo.params["lscale"] = 1e9
        geo.params["lpore"] = self.lpore
        self.add_synonymes(geo)
        self.add_curved_boundaries(geo)

    def add_synonymes(self, geo):
        # define porecurrent depending on molecule position
        porecurrent = "pore0"
        if porecurrent == self.where_is_molecule():
            porecurrent = "pore%d" % (self.nsections - 1,)
        dom = self.polygons[porecurrent]
        geo.params["lporecurrent"] = dom.b[1] - dom.a[1]

        poresolid = set(self.names)
        poresolidb = set([b for b in self.boundaries if any(
                         [b.startswith(name) for name in self.names])])
        pore = set([p for p in self.polygons if p.startswith("pore")])
        self.synonymes = dict(default_synonymes)
        self.synonymes.update(
            porecurrent = porecurrent,
            poresolid = poresolid,
            pore = pore,
            poresolidb = poresolidb,
        )
        geo.import_synonymes(self.synonymes)

    def add_curved_boundaries(self, geo):
        if self.dim == 2:
            if self.molecule:
                molec = curved.Circle(self.params.rMolecule,
                                      self.params.x0[::2])
                geo.curved = dict(moleculeb = molec.snap)
        elif self.dim == 3:
            raise NotImplementedError

    def to_gmsh(self, h=1.):
        # create mappings for nodes, edges
        self.gmsh_nodes = {}
        self.gmsh_edges = {}
        self.gmsh_cyl_surfs = {}
        self.gmsh_ball_surfs = {}

        lcs = self.set_length_scales(h)

        # create volumes
        for pname, p in self.balls.items() + self.polygons.items():
            gmsh.Comment("Creating %s subdomain." %pname)
            if isempty(p):
                gmsh.NoPhysicalVolume(pname)
                continue
            lc = lcs[pname]
            #print pname, lc
            self.Volume(p, pname, lc)

        # add physical surfaces
        for bname, bset in self.boundaries.items():
            self.PhysicalBoundary(bset, bname)

        gmsh.raw_code(["General.ExpertMode = 1;"])
        #gmsh.raw_code(["Mesh.Algorithm3D = 2;"])
        code = gmsh.get_code()
        meta = gmsh.get_meta()
        self.code = code
        self.meta = meta
        return code, meta

    def set_length_scales(self, h):
        lc = {p: h for p in self.polygons}
        lc.update({p: h*b.lc for p, b in self.balls.items() if not isempty(b)})
        lc["molecule"] = h*self.params.lcMolecule
        for i in range(self.nsections):
            lc["pore%d" % i] = h*self.params.lcCenter
        #lc.update(dict.fromkeys(self.names, h*self.params.lcCenter))
        return lc

    def Volume(self, p, pname, lc):
        dim = self.dim

        if dim == 2:
            surfs = [self.Edge(edge, lc) for edge in p.edges]
        elif dim == 3:
            surfs = self.Surfaces(p, lc)
            if hasattr(p, "holes"):
                for hole in p.holes:
                    surfs.extend(["-%s" % s for s in self.Surfaces(hole, lc)])
        else:
            raise NotImplementedError

        ll = FacetLoop[dim](surfs)
        vol = Entity[dim](ll)
        gmsh.PhysicalVolume(vol, pname, dim)

    def Surfaces(self, p, lc):
        if isinstance(p, Ball):
            return self.BallSurfaces(p, lc)
        elif isinstance(p, Polygon):
            return [s for e in p.edges for s in self.CylSurfaces(e, lc)]
        else:
            raise NotImplementedError

    def BallSurfaces(self, ball, lc):
        x0, r = tuple(ball.x0), ball.r
        if (x0, r) in self.gmsh_ball_surfs:
            return self.gmsh_ball_surfs[(x0, r)]

        surfs = add_ball(x0, r, lc)
        self.gmsh_ball_surfs[(x0, r)] = surfs
        return surfs

    def CylSurfaces(self, e, lc):
        if e in self.gmsh_cyl_surfs:
            return self.gmsh_cyl_surfs[e]

        # do nothing if edge is in center
        if e[0][0] == e[1][0] == 0.:
            return []

        # create Line
        line = self.Edge(e, lc)
        surfs = []

        # rotate in 4 steps
        for j in range(4):
            ex = gmsh.Extrude('Line{%s}' % line, rotation_axis=[0., 0., 1.],
                           point_on_axis=[0., 0., 0.], angle="Pi/2.0")
            line = ex + "[0]"
            surfs.append(ex + "[1]")

        # save, also for flipped edge
        self.gmsh_cyl_surfs[e] = surfs
        self.gmsh_cyl_surfs[e[::-1]] = ["-%s" % s for s in surfs]
        return surfs

    def Point(self, x, lc):
        if x in self.gmsh_nodes:
            return self.gmsh_nodes[x]

        #print x, lc
        if self.dim < 3:
            x1 = x + tuple(0. for i in range(3 - self.dim))
        else:
            x1 = [x[0], 0., x[1]]
        gmsh_x = gmsh.Point(x1, lc)
        self.gmsh_nodes[x] = gmsh_x
        return gmsh_x

    def Edge(self, e, lc=1.):
        # generalizes Line and Circle
        # if exists, return gmsh edge
        if e in self.gmsh_edges:
            return self.gmsh_edges[e]
        # otherwise, discriminate between Line and Circle
        points = [self.Point(v, lc) for v in e]

        if len(e) == 2:
            gmsh_e = gmsh.Line(*points)
        elif len(e) == 3:
            gmsh_e = gmsh.Circle(points)

        # save, also ffor flipped edge
        self.gmsh_edges[e] = gmsh_e
        self.gmsh_edges[e[::-1]] = "-" + gmsh_e
        return gmsh_e

    def PhysicalBoundary(self, bset, bname):
        gmsh.Comment("Creating %s boundary." %bname)
        dim = self.dim
        if not bset: # and not (dim == 3 and bname.startswith("molecule")):
            gmsh.NoPhysicalSurface(bname)
        elif dim == 2:
            boundary = [self.Edge(e) for e in bset]
            gmsh.PhysicalSurface(boundary, bname, dim)
        else:
            pass

class MultiPore(MultiPolygonPore, Pore):

    def __init__(self, polygons=None, **params):
        params = dict(default, **params)
        self.dim = params["dim"]
        self.geoname = params["geoname"] # name of folder where mesh is created
        MultiPolygonPore.__init__(self, **params)
        if polygons is not None:
            self.add_polygons(**polygons)

FacetLoop = {
    1: lambda x: x,
    2: gmsh.LineLoop,
    3: gmsh.SurfaceLoop
    }

Entity = {
    1: lambda x: gmsh.Line(*x),
    2: gmsh.PlaneSurface,
    3: gmsh.Volume
    }

def add_ball(m, r, lc):
    # add ball in 3D
    return gmsh.add_ball(m, r, lc, with_volume=False)[2]

def get_geo(poly, h=1., reconstruct=False, **params):
    p = Pore(poly, **params)

    if reconstruct:
        geo = maybe_reconstruct_geo(params=p.params)
        if geo is not None:
            p.build_polygons()
            p.build_boundaries()
            p.finalize_geo(geo)
            return geo

    geo = p.build(h=h)
    return geo

def maybe_reconstruct_geo(params=None):
    # return None if it does not work
    name = params["geoname"] if params is not None else "cylpore"
    try:
        geo = reconstructgeo(name=name, params=dict(params))
    except EnvironmentError as e:
        print e.message
        geo = None
    return geo

if __name__ == "__main__":
    from nanopores import user_params, showplots
    params = user_params(
        h = 1.,
        porename = "dna",
        H = 20.,
        R = 10.,
        cs = [1.7, -1.7],
        x0 = [0.,0.,0.],
    )

    dnapolygon = [[1, -5], [1, 5], [3, 5], [3, -5]]

    geo = get_geo(dnapolygon, **params)
    print geo
    print "params", geo.params

    geo.plot_subdomains()
    geo.plot_boundaries(interactive=True)

    # --- TEST
    p = MultiPore(dict(dna=dnapolygon), **params)
    p.build_polygons()
    p.build_boundaries()
    #print p.polygons.keys()
    #print p.boundaries.keys()

    p.protein.plot(".k", zorder=100)
    #p.polygons["bulkfluid_top"].plot()
    p.polygons["pore0"].plot()

    plot_edges(p.boundaries["memb"], color="blue")
    plot_edges(p.boundaries["lowerb"])
    plot_edges(p.boundaries["upperb"])
    plot_edges(p.boundaries["sideb"], color="yellow")
    plot_edges(p.boundaries["dnab"], color="red")
    showplots()
    # ---