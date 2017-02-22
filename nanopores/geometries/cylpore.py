# (c) 2017 Gregor Mitscha-Baude
import nanopores.py4gmsh as gmsh
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
    crosssections = (), # list of z coordinates
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
    crosssections = sorted(params.crosssections)

    # compute pore maximum and minimum z, length
    Z = [x[1] for x in poly]
    ztop = max(Z)
    zbot = min(Z)
    hpore = ztop - zbot
    assert all(zbot < z < ztop for z in crosssections)

    #Anchor Points on aHem for membran (index)
    ap1 = 18
    ap2 = 49
    apdiff=ap2-ap1

    #Anchor Points in aHem for CrossS (index)
    ac1 = 52
    ac2 = 68
    ac3 = 82
    ac4 = 0
    zcross = sorted([X_aHem[i][2] for i in [ac1, ac2, ac3, ac4]])
    params["lbtm"] = -zcross[0] + zcross[1]
    params["lctr"] = -zcross[1] + zcross[2]
    params["ltop"] = -zcross[2] + zcross[3]

    params["zporetop"] = zcross[3]
    params["zporebtm"] = zcross[0]
    params["ztop"] = params["zporetop"] + l3
    params["zbtm"] = params["zporebtm"] - l4

    r0=max([X_aHem[index][0] for index in [ac1, ac2, ac3, ac4]])+rMolecule

    X_Fluid_ext = numpy.array([[0.0, 0.0, l3],
                               [R, 0.0, l3],
                               [R, 0.0, X_aHem[ap1][2]],
                               [R, 0.0, X_aHem[ap2][2]],
                               [R, 0.0, -l0-l1-l4],
                               [0.0, 0.0, -l0-l1-l4]])
    X_Fluid_ctr = numpy.array([[0.0, 0.0, X_aHem[ac1][2]],
                               [0.0, 0.0, X_aHem[ac2][2]],
                               [0.0, 0.0, X_aHem[ac3][2]],
                               [0.0, 0.0, X_aHem[ac4][2]]])


    p_Fluid = [Point(x, lcOuter) for x in X_Fluid_ext]
    p_Fluid.extend([Point(y, lcCenter) for y in X_Fluid_ctr])
    p_aHem = [Point(x, lcCenter) for x in X_aHem]

    #Create Line Loops from the points sitting on the line
    Comment(' Connect all Fluid points ')
    e_Fluid = [Line(p_Fluid[k], p_Fluid[k+1]) for k in range(len(p_Fluid)-1)]
    e_Fluid.append(Line(p_Fluid[-1], p_Fluid[0]))

    Comment(' Connect all aHem points ')
    e_aHem = [Line(p_aHem[k], p_aHem[k+1]) for k in range(len(p_aHem)-1)]
    e_aHem.append(Line(p_aHem[-1], p_aHem[0]))

#    e_Membrane = [Line(p_aHem[ap1],p_Fluid[2]), Line(p_Fluid[3], p_aHem[ap2])]
    e_Membrane = [Line(Point(X_aHem[ap1],lcMembrane), Point(numpy.array([R,0.,X_aHem[ap1][2]]),lcMembrane)),\
                  Line(Point(numpy.array([R,0.,X_aHem[ap2][2]]),lcMembrane),Point(X_aHem[ap2],lcMembrane))]

    edges_to_rot = [e_Fluid[0:5], e_aHem, e_Membrane]

    geo_cs_str = "no crosssectional surface"
    if crosssections:
        Comment(' integrate crosssectional lines in fluid and check if molecule intersects lines')
        e_CrossS = [Line(p_aHem[ac1], p_Fluid[6]), Line(p_aHem[ac2], p_Fluid[7]), Line(p_aHem[ac3], p_Fluid[8]), Line(p_aHem[ac4], p_Fluid[9])]
        cs_pop_i = None
        # check if molecule is near pore
        if x0 is not None and (x0[0]**2 + x0[1]**2 <= r0**2):
            # check z coordinate of molecule
            if abs(x0[2] - X_aHem[ac4][2]) < rMolecule:
                geo_cs_str = "top crosssection"
                cs_pop_i = -1
            elif abs(x0[2] - X_aHem[ac3][2]) < rMolecule:
                geo_cs_str = "center top crosssection"
                cs_pop_i = 2
            elif abs(x0[2] - X_aHem[ac2][2]) < rMolecule:
                geo_cs_str = "center bottom crosssection"
                cs_pop_i = 1
            elif abs(x0[2] - X_aHem[ac1][2]) < rMolecule:
                geo_cs_str = "bottom crosssection"
                cs_pop_i = 0
        if cs_pop_i is not None:
            e_CrossS.pop(cs_pop_i)
            if cs_pop_i == 0:
                top_acdiff = len(X_aHem)-ap1
                bottom_end = ac2
            elif cs_pop_i == 1:
                top_acdiff = len(X_aHem)-ap1
                bottom_end = ac3
            elif cs_pop_i == 2:
                top_acdiff = ac2-ap1
                bottom_end = ac1
            elif cs_pop_i == -1:
                top_acdiff = ac3-ap1
                bottom_end = ac1
        edges_to_rot.append(e_CrossS)
        top_acdiff = len(X_aHem)-ap1
        bottom_end = ac1

if __name__ == "__main__":
    from alphahempoly import poly
    cross = [-0.5]
    get_geo(poly, h=1., crosssections=cross)