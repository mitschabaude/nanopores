# (c) 2017 Gregor Mitscha-Baude
import nanopores.py4gmsh as gmsh
from nanopores.tools.polygons import PolygonPore, plot_edges
from nanopores.tools.utilities import Params

default = dict(
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
    def compute_entities(self):
        # after build_...
        # create lists with common nodes, edges
        self.nodes = set([x for p in self.polygons.values() for x in p.nodes])
        self.edges = set([e for p in self.polygons.values() for e in p.edges])
        self.gmsh_nodes = [None for x in self.nodes]
        self.gmsh_edges = [None for x in self.edges]

    def entities_to_gmsh(self, dim=2):
        # create line loops for every polygon
        lineloops = {}
        for pname, p in self.polygons.items():
            ll = self.LineLoop(p.edges)
            lineloops[pname] = ll

        # if molecule, create circle
        protein = self.protein
        ll_protein = self.LineLoop(protein.edges)
        gmsh.PlaneSurface(ll_protein)



    def Point(self, x, lc):
        i = self.nodes.index(x)
        gmsh_x = self.gmsh_nodes[i]
        if gmsh_x is not None:
            return gmsh_x

        x = x + tuple(0. for i in range(3 - self.dim))
        gmsh_x = gmsh.Point(x, lc)
        self.gmsh_nodes[i] = gmsh_x
        return gmsh_x

    def Line(self, e, lc):
        # if exists, return edge
        i = self.edges.index(e)
        gmsh_e = self.gmsh_edges[i]
        if gmsh_e is not None:
            return gmsh_e

        x, y = e
        gmsh_x = self.Point(x, lc)
        gmsh_y = self.Point(y, lc)
        gmsh_e = gmsh.Line(gmsh_x, gmsh_y)
        self.gmsh_edges[i] = gmsh_e

        # also save flipped edge
        i1 = self.gmsh_edges.index((y, x))
        self.gmsh_edges[i1] = "-" + gmsh_e
        return gmsh_e

    def LineLoop(self, ll, lc):
        lines = [self.Line(edge, lc) for edge in ll]
        return gmsh.LineLoop(lines)



def get_geo(poly, h=1., **params):
    # get params
    params = Params(default, **params)
    p = Pore(poly, name="ahem", **params)
    p.build_polygons()
    p.build_boundaries()

    # --- TEST
    print p.polygons.keys()
    print p.boundaries.keys()

    p.protein.plot(".k", zorder=100)
    #p.polygons["bulkfluid_top"].plot()
    p.polygons["pore0"].plot()

    plot_edges(p.boundaries["memb"], color="blue")
    plot_edges(p.boundaries["lowerb"])
    plot_edges(p.boundaries["upperb"])
    plot_edges(p.boundaries["sideb"], color="yellow")
    plot_edges(p.boundaries["ahemb"], color="red")
    # ---






if __name__ == "__main__":
    from alphahempoly import poly
    from matplotlib import pyplot as plt
    cs = [-3, -6]
    geo = get_geo(poly, h=1., zmem=-5., cs=cs, proteincs=[-5.])
    plt.show()