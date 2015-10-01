"""
python script that generates mesh for Howorka geometry
"""

import numpy
from importlib import import_module
import nanopores.py4gmsh.basic
import nanopores.py4gmsh.extra
from nanopores.py4gmsh import *
from .params_geo import *

def get_geo(**params):
    """
    writes a 2d geo file for an axissymmetric geometry for Wei et al
    'Stochastic sensing ...'
    _________
    |        |
    |        |
    |   _____|
    |   \____|
    |        |
    |        |
    |________|

    """
    basic = reload(nanopores.py4gmsh.basic)
    extra = reload(nanopores.py4gmsh.extra)

    globals().update(params)
    Rx = R
    Ry = Rz

    # define additional geo variables
    # mesh generation should also work for no SAM layer
    if lsam < tolc or lsam is None:
        sam = None
    else:
        sam = True
    l0 =  lau +lsin +lsam

    angle2 = angle/2.0
    tan = numpy.tan(angle2*numpy.pi/180)
    sin = numpy.sin(angle2*numpy.pi/180)
    cos = numpy.cos(angle2*numpy.pi/180)
    r1 = r0 + l0*tan
    rsam = r0 + lsam/cos
    rsin = r0 + lsam/cos + rlau

    X_Fluid_up = numpy.array([
        [0, Ry, 0],
        [Rx, Ry, 0],
        [Rx, l0/2, 0],
    ])

    X_Fluid_low = numpy.array([
        [Rx, -l0/2, 0],
        [Rx, -Ry, 0],
        [0, -Ry, 0],
    ])

    X_Fluid_ctr = numpy.array([
        [0, -l0/2, 0],
        [0, -l0/6, 0],
        [0, +l0/6, 0],
        [0, +l0/2, 0],
    ])

    X_SAM_ctr = numpy.array([
        [r1, l0/2, 0],
        [(2*r1+r0)/3, l0/6, 0],
        [(r1+2*r0)/3, -l0/6, 0],
        [r0, -l0/2, 0],
    ])

    X_SiN = numpy.array([
        [Rx, -(l0/2-lsin), 0],
        [rsin +tan*lsin, -l0/2 +lsin, 0],
        [rsin, -l0/2, 0],
    ])

    p_SAM = [Point(x, lcCenter) for x in X_SAM_ctr]
    p_SiN = [Point(x, lcOuter) for x in X_SiN]
    p_Fluid_up = [Point(x, lcOuter) for x in X_Fluid_up]
    p_Fluid_low = [Point(x, lcOuter) for x in X_Fluid_low]
    p_Fluid_ctr = [Point(x, lcCenter) for x in X_Fluid_ctr]

    if sam is None:
        p_Au = []
    else:
        X_Au = numpy.array([
            [Rx, -lsam +l0/2, 0],
            [rsam -tan*(lsam-l0), -lsam + l0/2, 0],
            [rsam, -l0/2, 0],
        ])
        p_Au = [Point(X_Au[0], lcOuter)]
        p_Au.extend([Point(x, lcCenter) for x in X_Au[1:]])

    # Group all fluid points into list
    p_Fluid = p_Fluid_up + p_SAM
    if sam is not None:
        p_Fluid.append(p_Au[-1])
    p_Fluid.append(p_SiN[-1])
    p_Fluid.extend(p_Fluid_low + p_Fluid_ctr)

    geo_cs_str = "no crosssectional surface"
    cs_pop_i = None
    insert_ind = 0

    Comment(' integrate crosssectional lines in fluid surface ')
    z_CrossS = [X_Fluid_ctr[k][1] for k in reversed(range(len(X_Fluid_ctr)))]
    e_CrossS = [Line(p_Fluid_ctr[k], p_SAM[len(p_SAM)-1-k]) for k in reversed(range(len(p_Fluid_ctr)))]

    if x0 is not None:
        X_Molecule = numpy.array([[0, x0[2] -rMolecule, 0],
                                  [0, x0[2], 0],
                                  [0, x0[2] +rMolecule, 0]])
        p_Molecule = [Point(x, lcMolecule) for x in X_Molecule]

        geo_cs_list = ["top", "center top", "center bottom", "bottom"]
        # Determine position of molecule for correct geometry
        for i in range(len(z_CrossS)):
            if abs(x0[2] - z_CrossS[i]) < rMolecule:
                cs_pop_i = i
            if x0[2] < z_CrossS[i] - rMolecule:
                insert_ind = i+1

    if cs_pop_i is not None:
        geo_cs_str = geo_cs_list[cs_pop_i]
        e_CrossS.pop(cs_pop_i)
        p_Fluid.pop(len(p_Fluid)-1-cs_pop_i)

    if x0 is not None:
        p_Fluid.insert(len(p_Fluid)-insert_ind, p_Molecule[0])
        p_Fluid.insert(len(p_Fluid)-insert_ind, p_Molecule[2])
        c_Molecule = Circle(p_Molecule)
        e_Molecule = Line(p_Molecule[2], p_Molecule[0])
        ll_Molecule = LineLoop([c_Molecule, e_Molecule])
        s_Molecule = PlaneSurface(ll_Molecule)

    # Create Line Loops from the points sitting on the line
    Comment(' Connect all Fluid points ')
    e_Fluid = [Line(p_Fluid[k], p_Fluid[k+1]) for k in range(len(p_Fluid)-1)]
    e_Fluid.append(Line(p_Fluid[-1], p_Fluid[0]))
    if x0 is not None:
        e_Fluid.pop(len(e_Fluid) - insert_ind -2)
        e_Fluid.insert(len(e_Fluid) - insert_ind -1, c_Molecule)
    ll_Fluid = LineLoop(e_Fluid)

    Comment(' Connect all outer Membrane points ')
    slicem = slice(len(p_Fluid_up) -1, p_Fluid.index(p_SiN[-1]) +1)
    e_Membrane = e_Fluid[slicem]
    if sam is not None:
        p_list_mem = [p_Fluid_low[0], p_SiN[0], p_Au[0], p_Fluid_up[-1]]
    else:
        p_list_mem = [p_Fluid_low[0], p_SiN[0], p_Fluid_up[-1]]
    e_Membrane.extend([Line(p_list_mem[k], p_list_mem[k+1])  \
                       for k in range(len(p_list_mem)-1)])
    ll_Membrane = LineLoop(e_Membrane)

    s_Fluid = PlaneSurface(ll_Fluid)
    s_Membrane = PlaneSurface(ll_Membrane)

    Comment(' Integrate crossections into fluid surface')
    raw_code(['Line{%s} In Surface{%s};'  \
              %(e_CrossS[k], s_Fluid) for k in range(len(e_CrossS))])

    Comment(' Integrate Membrane material interfaces into membrane surface ')
    if sam is not None:
        e_Au = [Line(p_Au[k], p_Au[k+1]) for k in range(len(p_Au)-1)]
        raw_code(['Line{%s} In Surface{%s};'  \
                  %(e_Au[k], s_Membrane) for k in range(len(e_Au))])

    e_SiN = [Line(p_SiN[k], p_SiN[k+1]) for k in range(len(p_SiN)-1)]
    raw_code(['Line{%s} In Surface{%s};'  \
              %(e_SiN[k], s_Membrane) for k in range(len(e_SiN))])


    # Define Fields for varying mesh size
    field_list = []
    if membraneblayer == True:
        blayer = BoundaryLayer(edges_list=e_Membrane[1:5], hfar=lcOuter,
                               hwall_n=lcCenter*0.2, hwall_t=lcCenter*0.5,
                               thickness=1, ratio=2)
        field_list.append(blayer)

    if boxfields:
        box1 = BoxField(lcCenter, lcOuter, -1e-15, r1*1.3, \
                        -l0/2 -r0, l0/2 +r0)
        field_list.append(box1)

    if lowerblayer:
        n_edge_l = 9 if not sam else 10
        blayer_l = BoundaryLayer(
            edges_list=e_Fluid[n_edge_l:n_edge_l+1],
            hfar=lcOuter, hwall_n=lcCenter*0.5, hwall_t=lcCenter*1.0,
            thickness=1, ratio=2)
        field_list.append(blayer_l)

    if x0 is not None and moleculeblayer:
        blayer_mol = BoundaryLayer(
            edges_list=[c_Molecule,], hfar=lcOuter,
            hwall_n=lcCenter*0.1, hwall_t=lcCenter*0.5,
            thickness=1, ratio=2)
        field_list.append(blayer_mol)

    raw_code(['bfield = newf;'])
    raw_code(['Field[bfield] = Min;'])
    raw_code(['Field[bfield].FieldsList = {%s};' %(','.join(field_list))])
    raw_code(['Background Field = bfield;'])

    # to disable question dialogs
    raw_code(['General.ExpertMode = 1;'])
    # Don't extend the elements sizes from the boundary inside the domain
    raw_code(['Mesh.CharacteristicLengthExtendFromBoundary = 0;'])
    # 2D algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
    # only delaunay is compatible with attractor fields
    raw_code(['Mesh.Algorithm = 5;'])

    geo_dict = {"gmsh mesh generating sript": __name__,
                "xMolecule": x0,
                "Number of crosssections": len(e_CrossS),
                "Total number of crossections": 4,
                "molecule crosses": geo_cs_str,
                "popped crossection index": cs_pop_i,
                "cs_pop_i": cs_pop_i,
                "Typical length scale": lc,
                "geo_code": get_code(),
            }
    return geo_dict


# -----
if __name__ == '__main__':
    print(get_geo())
    print('\n - This is the sample code for the geo file')
