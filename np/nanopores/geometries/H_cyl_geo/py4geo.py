import numpy, math
from .params_geo import *
import nanopores.py4gmsh.basic
import nanopores.py4gmsh.extra
from nanopores.py4gmsh import *
from warnings import warn

def get_geo(x0 = None, crosssections = True, **params):
    """
    writes a 3d geo file for an extruded axissymmetric geometry
    for Howorka 'Self-Assembled DNA that spans lipid bilayers'
    _________
    |        |
    |  _     |
    | | |____|
    | | |____|
    | |_|    |
    |        |
    |________|   *rotated around y-axis *

    """
    reload(nanopores.py4gmsh.basic)
    reload(nanopores.py4gmsh.extra)
    globals().update(params)

    X_Fluid_ext = numpy.array([[0.0, 0.0, Rz],
                               [R, 0.0, Rz],
                               [R, 0.0, l1/2],
                               [R, 0.0, -l1/2],
                               [R, 0.0, -Rz],
                               [0.0, 0.0, -Rz]])
    X_Fluid_ctr = numpy.array([[0.0, 0.0, -l0/2],
                               [0.0, 0.0, -l1/2],
                               [0.0, 0.0, l1/2],
                               [0.0, 0.0, l0/2]])
    X_DNA = numpy.array([[r0, 0.0, l0/2],
                         [r1, 0.0, l0/2],
                         [r1, 0.0, l1/2],
                         [r1, 0.0, -l1/2],
                         [r1, 0.0, -l0/2],
                         [r0, 0.0, -l0/2],
                         [r0, 0.0, -l1/2],
                         [r0, 0.0, l1/2]])

    p_Fluid = [Point(x, lcOuter) for x in X_Fluid_ext]
    p_Fluid.extend([Point(y, lcCenter) for y in X_Fluid_ctr])
    p_DNA = [Point(x, lcCenter) for x in X_DNA]

    #Create Line Loops from the points sitting on the line
    Comment(' Connect all Fluid points ')
    e_Fluid = [Line(p_Fluid[k], p_Fluid[k+1]) for k in range(len(p_Fluid)-1)]
    e_Fluid.append(Line(p_Fluid[-1], p_Fluid[0]))

    Comment(' Connect all DNA points ')
    e_DNA = [Line(p_DNA[k], p_DNA[k+1]) for k in range(len(p_DNA)-1)]
    e_DNA.append(Line(p_DNA[-1], p_DNA[0]))

    e_Membrane = [Line(p_DNA[2], p_Fluid[2]), Line(p_Fluid[3], p_DNA[3])]

    edges_to_rot = [e_Fluid[0:5], e_DNA, e_Membrane]

    geo_cs_str = "no crosssectional surface"
    if crosssections:
        Comment(' integrate crosssectional lines in fluid and check if molecule intersects lines')
        e_CrossS = [Line(p_DNA[k], p_Fluid[k+1]) for k in range(5,len(p_DNA))]
        e_CrossS.append(Line(p_DNA[0], p_Fluid[-1]))
        cs_pop_i = None
        # check if molecule is near pore
        if x0 is not None and (x0[0]**2 + x0[1]**2 <= r0**2):
            # check z coordinate of molecule
            if abs(x0[2] - l0/2) < rMolecule:
                geo_cs_str = "top crosssection"
                cs_pop_i = -1
            elif abs(x0[2] - l1/2) < rMolecule:
                geo_cs_str = "center top crosssection"
                cs_pop_i = 2
            elif abs(x0[2] + l1/2) < rMolecule:
                geo_cs_str = "center bottom crosssection"
                cs_pop_i = 1
            elif abs(x0[2] + l0/2) < rMolecule:
                geo_cs_str = "bottom crosssection"
                cs_pop_i = 0
        if cs_pop_i is not None:
            e_CrossS.pop(cs_pop_i)
        edges_to_rot.append(e_CrossS)

    rot_axis = [0.0, 0.0, 1.0]
    point_on_rot_axis = [0.0, 0.0, 0.0]
    # Extrude all edge 4 times Pi/2
    surfs = []
    angle = 'Pi/2'
    n_e = len(edges_to_rot)
    n_e_i = [len(edges_to_rot[i]) for i in range(n_e)]
    for i in range(n_e):
        surfs_i = []
        Comment('Extrude in 4 steps around z-axis.')
        previous = edges_to_rot[i]
        for j in range(4):
            Comment('Step %s' % (j+1))
            for k in range(len(previous)):
                name = Extrude('Line{%s}' % previous[k],
                               rotation_axis=rot_axis,
                               point_on_axis=point_on_rot_axis,
                               angle=angle
                           )
                surfs_i.append(name + '[1]')
                previous[k] = name + '[0]'
        surfs.append(surfs_i)

    surfs_Fluid = surfs[0][:]  # [:] is important for a shallow copy (-> del nextline)
    del surfs_Fluid[2::n_e_i[0]]  # deletes outer membrane boundary
    surfs_Fluid_DNA = surfs[1][:]
    del surfs_Fluid_DNA[2::n_e_i[1]]  # deletes membrane
    [surfs_Fluid.append(s) for s in surfs_Fluid_DNA]
    [surfs_Fluid.append(s) for s in surfs[2]]
    sl_Fluid = SurfaceLoop(surfs_Fluid)

    sl_DNA = SurfaceLoop(surfs[1])
    vol_DNA = Volume(sl_DNA)

    surfs_Membrane = surfs[0][2::n_e_i[0]]
    [surfs_Membrane.append(s) for s in surfs[1][2::n_e_i[1]]]
    [surfs_Membrane.append(s) for s in surfs[2]]
    sl_Membrane = SurfaceLoop(surfs_Membrane)
    vol_Membrane = Volume(sl_Membrane)

    if x0 is None:
        vol_Fluid = Volume(sl_Fluid)
    else:
        Comment('Add molecule ball')
        Molecule = add_ball(numpy.asarray(x0), rMolecule, lcMolecule,
                            with_volume=True, holes=None, label=None
                        )
        sl_Fluid_Molecule = Array([sl_Fluid] + [Molecule[1]])
        vol_Fluid = Volume(sl_Fluid_Molecule)
        # Molecule[0]->Volume,  Molecule[1]->surface loop,  Molecule[2]->surfs
        vol_Molecule = Molecule[0]

    if crosssections:
        surfs_CrossS = surfs[3]
        raw_code(['Surface{%s} In Volume{%s};'  \
                  %(surfs_CrossS[k], vol_Fluid) for k in range(len(surfs_CrossS))])

    if membraneblayer == True:
        warn('Currently no membrane boundary layers implemented for this geometry')
    membraneblayer_list = []

    if moleculeblayer and x0 is not None:
        moleculeblayer_list = Molecule[2]
    else:
        moleculeblayer_list = []

    blayer_list = moleculeblayer_list + membraneblayer_list
    if blayer_list:
        blayer = BoundaryLayer(
            edges_list=None, faces_list=blayer_list,
            hfar=lcOuter, hwall_n=lcCenter*0.1, hwall_t=lcCenter*0.5,
            thickness=1, ratio=2)
        field_list = [blayer,]

        raw_code(['bfield = newf;'])
        raw_code(['Field[bfield] = Min;'])
        raw_code(['Field[bfield].FieldsList = {%s};' %(','.join(field_list))])
        # Uncomment for mesh size field
        raw_code(['Background Field = bfield;'])

    # to disable question dialogs
    raw_code(['General.ExpertMode = 1;'])

    # Meshing Algorithm: 2= ?, 5 = frontal (netgen)
    #raw_code(['Mesh.Algorithm3D = 5;'])

    geo_dict = {"gmsh mesh generating sript": __name__,
                "xMolecule": x0,
                "crosssections": crosssections,
                "Number of crosssections": len(e_CrossS),
                "Total number of crossections": 4,
                "molecule crosses": geo_cs_str,
                "popped crossection index": cs_pop_i,
                "cs_pop_i": cs_pop_i,
                "Typical length scale on DNA": lcCenter,
                "geo_code": get_code(),
            }
    return geo_dict


# -----
if __name__ == '__main__':
    print(get_geo())
    print('\n - This is the sample code for the geo file')
