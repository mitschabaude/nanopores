# (c) 2017 Gregor Mitscha-Baude
import nanopores.py4gmsh as gmsh
import nanopores.geometries.curved as curved
from nanopores.geo2xml import geofile2geo, reconstructgeo
from nanopores.tools.polygons import (Ball, Polygon,
                                      PolygonPore, MultiPolygonPore,
                                      plot_edges, isempty)
from nanopores.tools.utilities import Log

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
    cs = None, # crosssections; list of z coordinates
    poreregion = False, # whether to include fluid above pore as subdomain
)

default_synonymes = dict(
    #subdomains
    bulkfluid = {"bulkfluid_top", "bulkfluid_bottom"},
    fluid = {"bulkfluid", "pore", "nearpore"},
    nearpore = {"poreregion_top", "poreregion_bottom"},
    poreregion = {"pore", "nearpore"},
    solid = {"membrane", "poresolid", "molecules"},
    ions = "fluid",
    #boundaries
    noslip = {"poresolidb", "memb", "moleculesb"},
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

    def build(self, h=1., subs=None, reconstruct=False):
        with Log("computing polygons..."):
            self.build_polygons()
            self.build_boundaries()
        geo = self.build_geometry(h, subs, reconstruct)
        return geo

    def build_geometry(self, h=1., subs=None, reconstruct=False):
        self.add_synonymes()
        self.choose_domains(subs)
        self.geo = None
        if reconstruct:
            self.geo = maybe_reconstruct_geo(params=self.params)
        if self.geo is None:
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
        geo.import_synonymes(self.synonymes)
        self.add_curved_boundaries(geo)

    def add_synonymes(self):
        # define porecurrent depending on molecule position
        porecurrent = "pore0"
        if porecurrent == self.where_is_molecule():
            porecurrent = "pore%d" % (self.nsections - 1,)
        dom = self.polygons[porecurrent]
        self.params["lporecurrent"] = dom.b[1] - dom.a[1]

        poresolid = set(self.names)
        poresolidb = set([b for b in self.boundaries if any(
                         [b.startswith(name) for name in self.names])])
        pore = set(["pore%d" % i for i in range(self.nsections)])
        molecules = set(self.balls.keys())
        moleculesb = set(name + "b" for name in self.balls)
        self.synonymes = dict(default_synonymes)
        self.synonymes.update(
            porecurrent = porecurrent,
            poresolid = poresolid,
            pore = pore,
            poresolidb = poresolidb,
            molecules = molecules,
            moleculesb = moleculesb,
        )
    
    def unpack_synonymes(self, syns):
        if isinstance(syns, str):
            syns = {syns}
        domains = set()
        for syn in syns:
            if syn in self.domains:
                domains.add(syn)
            elif syn in self.synonymes:
                domains |= self.unpack_synonymes(self.synonymes[syn])
        return domains
    
    def choose_domains(self, subs):
        if subs is None:
            return
        domains = self.unpack_synonymes(subs)
        for dom in self.domains:
            if not dom in domains:
                self.domains.pop(dom)
        
    def add_curved_boundaries(self, geo):
        geo.curved = {}
        for bname, ball in self.balls.items():
            if not isempty(ball):
                name = bname + "b"
                if self.dim == 2:
                    subdomain = curved.Circle(ball.r, ball.x0)
                elif self.dim == 3:
                    subdomain = curved.Sphere(ball.r, ball.x0)
                geo.curved[name] = subdomain.snap     

    def to_gmsh(self, h=1.):
        # create mappings for nodes, edges
        self.gmsh_nodes = {}
        self.gmsh_edges = {}
        self.gmsh_cyl_surfs = {}
        self.gmsh_ball_surfs = {}

        lcs = self.set_length_scales(h)

        # create volumes
        for pname, p in self.domains.items():
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

    def BallSurfaces(self, ball, lc=1.):
        x0, r = tuple(ball.x0), ball.r
        if (x0, r) in self.gmsh_ball_surfs:
            return self.gmsh_ball_surfs[(x0, r)]

        surfs = add_ball(x0, r, lc)
        self.gmsh_ball_surfs[(x0, r)] = surfs
        return surfs

    def CylSurfaces(self, e, lc=1.):
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
    
    def get_boundary(self, e):
        dim = self.dim
        if dim == 2:
            if e in self.gmsh_edges:
                return [self.gmsh_edges[e]]
        elif dim == 3:
            if isinstance(e, Ball):
                bdict = self.gmsh_ball_surfs
                e = tuple(e.x0), e.r
            else:
                bdict = self.gmsh_cyl_surfs
            if e in bdict:
                return bdict[e]
        return []            

    def PhysicalBoundary(self, bset, bname):
        # should not create any new stuff
        gmsh.Comment("Creating %s boundary." %bname)
        boundary = [b for e in bset for b in self.get_boundary(e)]
        if len(boundary) > 0:
            gmsh.PhysicalSurface(boundary, bname, self.dim)
        else:
            gmsh.NoPhysicalSurface(bname)

class MultiPore(MultiPolygonPore, Pore):

    def __init__(self, polygons=None, balls=None, **params):
        params = dict(default, **params)
        self.dim = params["dim"]
        self.geoname = params["geoname"] # name of folder where mesh is created
        MultiPolygonPore.__init__(self, **params)
        if polygons is not None:
            self.add_polygons(**polygons)
        if balls is not None:
            self.add_balls(**balls)

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

def get_geo(poly, h=1., reconstruct=False, subs=None, **params):
    pore = Pore(poly, **params)
    geo = pore.build(h=h, reconstruct=reconstruct, subs=subs)
    return geo

def get_geo_multi(polys, balls=None, h=1., reconstruct=False, subs=None, **params):
    pore = MultiPore(polys, balls, **params)
    geo = pore.build(h=h, reconstruct=reconstruct, subs=subs)
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
        x0 = [0.,0.,0.],
        dim = 2,
        subs = None,
    )

    dnapolygon = [[1, -5], [1, 5], [3, 5], [3, -5]]
    
        # --- TEST
    p = MultiPore(dict(dna=dnapolygon), **params)
    p.build_polygons()
    p.build_boundaries()
    #print p.polygons.keys()
    #print p.boundaries.keys()

    p.protein.plot(".k", zorder=100)
    #p.polygons["bulkfluid_top"].plot()
    p.polygons["pore0"].plot()
    p.polygons["bulkfluid_top"].plot("--k")
    p.polygons["bulkfluid_bottom"].plot("--k")

    plot_edges(p.boundaries["memb"], color="blue")
    plot_edges(p.boundaries["lowerb"], color="green")
    plot_edges(p.boundaries["upperb"], color="green")
    plot_edges(p.boundaries["sideb"], color="yellow")
    plot_edges(p.boundaries["dnab"], color="red")
    showplots()
    # ---

    geo = get_geo_multi(dict(dna=dnapolygon), **params)
    print geo
    print "params", geo.params

    geo.plot_subdomains()
    geo.plot_boundaries(interactive=True)
