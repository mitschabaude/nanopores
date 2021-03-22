"""
python script that generates mesh for Howorka geometry
"""

from .params_geo import *
import numpy
from importlib import import_module
import nanopores.py4gmsh.basic
import nanopores.py4gmsh.extra
from nanopores.py4gmsh import *
from warnings import warn
import importlib

def get_geo(**params):
    """
    writes a 2d geo file for an axissymmetric geometry for Howorka
    'Self-Assembled DNA that spans lipid bilayers'
    _________
    |        |
    |  _     |
    | | |____|
    | | |____|
    | |_|    |
    |        |
    |________|

    """
    importlib.reload(nanopores.py4gmsh.basic)
    importlib.reload(nanopores.py4gmsh.extra)

    # characteristic length / mesh size h updated from params_geo
    globals().update(params)

    X_Fluid_ext = numpy.array([[0, Ry, 0],
                               [Rx, Ry, 0],
                               [Rx, l1/2, 0],
                               [Rx, -l1/2, 0],
                               [Rx, -Ry, 0],
                               [0, -Ry, 0], ])
    X_Fluid_ctr = numpy.array([[0, -l0/2, 0],
                               [0, -l1/2, 0],
                               [0, l1/2, 0],
                               [0, l0/2, 0], ])

    X_DNA = numpy.array([[r0, l0/2, 0],
                         [r1, l0/2, 0],
                         [r1, l1/2, 0],
                         [r1, -l1/2, 0],
                         [r1, -l0/2, 0],
                         [r0, -l0/2, 0],
                         [r0, -l1/2, 0],
                         [r0, l1/2, 0]])

    p_Fluid = [Point(x, lcOuter) for x in X_Fluid_ext]
    p_Fluid.extend([Point(x, lcCenter) for x in X_Fluid_ctr])
    p_DNA = [Point(x, lcCenter) for x in X_DNA]

    geo_cs_str = "no crosssectional surface"
    cs_pop_i = None
    insert_ind = 0

    Comment(' integrate crosssectional lines in fluid surface ')
    z_CrossS = [X_Fluid_ctr[k][1] for k in reversed(list(range(len(X_Fluid_ctr))))]
    e_CrossS = [Line(p_Fluid[9], p_DNA[0])]
    e_CrossS.extend([Line(p_Fluid[k], p_DNA[k-1]) for k in range(8,5,-1)])

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

    Comment(' Connect all DNA points ')
    e_DNA = [Line(p_DNA[k], p_DNA[k+1]) for k in range(len(p_DNA)-1)]
    e_DNA.append(Line(p_DNA[-1], p_DNA[0]))
    ll_DNA = LineLoop(e_DNA)

    Comment(' Membrane line loop requires extra work ')
    e_Membrane = [Line(p_DNA[2], p_Fluid[2]), Line(p_Fluid[3], p_DNA[3])]
    ll_Membrane = LineLoop([e_Membrane[0], e_Fluid[2], e_Membrane[1], '-%s' %e_DNA[2]])

    # Workaround PlaneSurface accepts only one LineLoop
    s_Fluid = PlaneSurface('%s' %','.join([ll_Fluid, ll_DNA, ll_Membrane]))
    s_DNA = PlaneSurface(ll_DNA)
    s_Membrane = PlaneSurface(ll_Membrane)

    raw_code(['Line{%s} In Surface{%s};'  \
              %(e_CrossS[k], s_Fluid) for k in range(len(e_CrossS))])

    # Define Fields for varying mesh size
    field_list = []
    if boxfields:
        box1 = BoxField(lcOuter, lcOuter, r0*0.4, r0*1.3, \
                    -l0*0.5-r0*0.4, l0*0.5+r0*0.4)
        box2 = BoxField(lcOuter, lcOuter, r1-r0*0.3, r1+r0*0.6,  \
                    -l0*0.5-r0*0.4, l0*0.5+r0*0.4)
        box3 = BoxField(lcCenter, lcOuter, -1e-15, r1*1.2, \
                    -l0*0.55, l0*0.55)
        field_list.extend([box1, box2, box3])

    if x0 is not None and boxfields:
        box4 = BoxField(lcMolecule, lcOuter, -1e-15, r0,
                        x0[2]-1.8*rMolecule, x0[2]+1.8*rMolecule)
        field_list.append(box4)

    if x0 is not None and moleculeblayer:
        blayer_mol = BoundaryLayer(
            edges_list=[c_Molecule,], hfar=lcOuter,
            hwall_n=lcCenter*0.1, hwall_t=lcCenter*0.5,
            thickness=1, ratio=2)
        field_list.append(blayer_mol)

    if membraneblayer:
        warn('Membrane boundary layers not defined for this geometry')

    if field_list:
        raw_code(['bfield = newf;'])
        raw_code(['Field[bfield] = Min;'])
        raw_code(['Field[bfield].FieldsList = {%s};' %', '.join(field_list)])
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
    print((get_geo()))
    print('\n - This is the sample code for the geo file')
