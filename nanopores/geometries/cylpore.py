# (c) 2017 Gregor Mitscha-Baude
import nanopores.py4gmsh as gmsh
from nanopores.tools.polygons import PorePolygon
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

def get_geo(poly, h=1., **params):
    # get params
    params = Params(default, **params)
    dim = params.dim
    R = params.R
    H = params.H
    x0 = params.x0
    rMolecule = params.rMolecule
    lcMolecule = params.lcMolecule
    lcCenter = params.lcCenter
    hmem = params.hmem
    zmem = params.zmem
    cs = sorted(params.cs)

    # compute pore maximum and minimum z, length
    Z = [x[1] for x in poly]
    ztop = max(Z)
    zbot = min(Z)

    poly = PorePolygon(poly)
    pmem = p.get_membrane(hmem, zmem, R)
    pbot, pctr, ptop = p.get_poresections(cs)

    p.plot()
    pmem.plot()
    top, bottom = p.get_bulkfluids(R, H, p, pmem, pctr, pbot)
    top.plot()
    bottom.plot()
    plt.show()



if __name__ == "__main__":
    from alphahempoly import poly
    cross = [-3, -6]
    get_geo(poly, h=1., crosssections=cross)