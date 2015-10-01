"""
python script that generates mesh for Wei et al geometry
"""

import numpy
import nanopores.py4gmsh.basic
import nanopores.py4gmsh.extra
from nanopores.py4gmsh import *
from params_geo import *
from warnings import warn

def get_geo(x0 = None, **params):
    """
    writes a 3d geo file for an axissymmetric geometry for Wei et al
    'Stochastic sensing ...'
    _________
    |        |
    |        |
    |   _____|
    |   \____|
    |        |
    |        |
    |________|   *rotated around z-axis

    """

    reload(nanopores.py4gmsh.basic)
    reload(nanopores.py4gmsh.extra)
    globals().update(params)

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
    Rx = R

    X_Fluid_up = numpy.array([
        [0, 0, Rz],
        [Rx, 0, Rz],
        [Rx, 0, l0/2],
    ])

    X_Fluid_low = numpy.array([
        [Rx, 0, -l0/2],
        [Rx, 0, -Rz],
        [0, 0, -Rz],
    ])

    X_Fluid_ctr = numpy.array([
        [0, 0, -l0/2],
        [0, 0, -l0/6],
        [0, 0, +l0/6],
        [0, 0, +l0/2],
    ])

    X_SAM_ctr = numpy.array([
        [r1, 0, l0/2],
        [(2*r1+r0)/3, 0, l0/6],
        [(r1+2*r0)/3, 0, -l0/6],
        [r0, 0, -l0/2],
    ])

    X_SiN = numpy.array([
        [Rx, 0, -(l0/2-lsin)],
        [rsin +tan*lsin, 0, -l0/2 +lsin],
        [rsin, 0, -l0/2],
    ])

    p_SAM = [Point(x, lcCenter) for x in X_SAM_ctr[:-1]]
    p_SAM.append(Point(X_SAM_ctr[-1],lcCenter/5.0))
    p_SiN = [Point(x, lcOuter) for x in X_SiN]
    p_Fluid_up = [Point(x, lcOuter) for x in X_Fluid_up]
    p_Fluid_low = [Point(x, lcOuter) for x in X_Fluid_low]
    p_Fluid_ctr = [Point(x, lcCenter) for x in X_Fluid_ctr]

    if sam is None:
        p_Au = []
    else:
        X_Au = numpy.array([
            [Rx, 0, -lsam +l0/2],
            [rsam -tan*(lsam-l0), 0, -lsam + l0/2],
            [rsam, 0, -l0/2],
        ])
        p_Au = [Point(X_Au[0], lcOuter)]
        p_Au.extend([Point(x, lcCenter) for x in X_Au[1:]])

    # Group all fluid points into list
    p_Fluid = p_Fluid_up + p_SAM
    if sam is not None:
        p_Fluid.append(p_Au[-1])
    p_Fluid.append(p_SiN[-1])
    p_Fluid.extend(p_Fluid_low)

    geo_cs_str = "no crosssectional surface"
    cs_pop_i = None
    insert_ind = 0

    Comment(' integrate crosssectional lines in fluid surface ')
    z_CrossS = [X_Fluid_ctr[k][2] for k in reversed(range(len(X_Fluid_ctr)))]
    e_CrossS = [Line(p_Fluid_ctr[k], p_SAM[len(p_SAM)-1-k]) for k in reversed(range(len(p_Fluid_ctr)))]

    if x0 is not None and (x0[0]**2 + x0[1]**2 <= r1**2):
        geo_cs_list = ["top", "center top", "center bottom", "bottom"]
        # Determine position of molecule for correct geometry
        for i in range(len(z_CrossS)):
            if abs(x0[2] - z_CrossS[i]) < rMolecule:
                cs_pop_i = i

    if cs_pop_i is not None:
        geo_cs_str = geo_cs_list[cs_pop_i]
        e_CrossS.pop(cs_pop_i)

    # Create Line Loops from the points sitting on the line
    Comment(' Connect all outer Fluid points ')
    e_Fluid = [Line(p_Fluid[k], p_Fluid[k+1]) for k in range(len(p_Fluid)-1)]

    Comment(' Connect outer Membrane points ')
    if sam is not None:
        p_list_mem = [p_Fluid_low[0], p_SiN[0], p_Au[0], p_Fluid_up[-1]]
    else:
        p_list_mem = [p_Fluid_low[0], p_SiN[0], p_Fluid_up[-1]]
    e_Membrane_ext = [Line(p_list_mem[k], p_list_mem[k+1])  \
                       for k in range(len(p_list_mem)-1)]

    e_Membrane_in = []
    Comment(' Integrate Membrane material interfaces into membrane')
    if sam is not None:
        e_Au = [Line(p_Au[k], p_Au[k+1]) for k in range(len(p_Au)-1)]
        e_Membrane_in.extend(e_Au)
    e_SiN = [Line(p_SiN[k], p_SiN[k+1]) for k in range(len(p_SiN)-1)]
    e_Membrane_in.extend(e_SiN)

    edges_to_rot = [e_Fluid, e_Membrane_ext, e_Membrane_in, e_CrossS]

    rot_axis = [0.0, 0.0, 1.0]
    point_on_rot_axis = [0.0, 0.0, 0.0]
    surfs = []
    n_rot = 4
    rot_angle = '2*Pi/%f' %n_rot
    n_e = len(edges_to_rot)
    n_e_i = [len(edges_to_rot[i]) for i in range(n_e)]
    for i in range(n_e):
        surfs_i = []
        Comment('Extrude in 4 steps by rot_angle %s around z-axis.' %rot_angle)
        previous = edges_to_rot[i]
        for j in range(n_rot):
            Comment('Step %s' % (j+1))
            for k in range(len(previous)):
                name = Extrude('Line{%s}' % previous[k],
                               rotation_axis=rot_axis,
                               point_on_axis=point_on_rot_axis,
                               angle=rot_angle
                           )
                surfs_i.append(name + '[1]')
                previous[k] = name + '[0]'
        surfs.append(surfs_i)

    surfs_Fluid = surfs[0][:]
    sl_Fluid = SurfaceLoop(surfs_Fluid)

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

    Comment(' Integrate crosssections into fluid')
    surfs_CrossS = surfs[-1]
    raw_code(['Surface{%s} In Volume{%s};' %(surfs_CrossS[k], vol_Fluid)  \
              for k in range(len(surfs_CrossS))])

    Comment(' Add membrane')
    surfs_Membrane = surfs[1][:]
    n_Membrane_Fluid = len(p_Fluid_up+p_SAM)-1
    for k in range(n_Membrane_Fluid):
        surfs_Membrane.extend(surfs[0][len(p_Fluid_up)+k-1::n_e_i[0]])
    sl_Membrane = SurfaceLoop(surfs_Membrane)
    vol_Membrane = Volume(sl_Membrane)

    Comment(' Integrate membrane interfaces into membrane')
    surfs_MembraneIn = surfs[-2]
    raw_code(['Surface{%s} In Volume{%s};' %(surfs_MembraneIn[k], vol_Membrane)  \
              for k in range(len(surfs_MembraneIn))])

    if moleculeblayer and x0 is not None:
        moleculeblayer_list = Molecule[2]
    else:
        moleculeblayer_list = []

    if membraneblayer:
        n_bl_start = len(surfs[1])+1*n_rot
        #numbers of surfaces with a boundary layer: 4, or if sam 5
        num_faces = (4 if not sam else 5)
        membraneblayer_list = surfs_Membrane[n_bl_start:n_bl_start+num_faces*n_rot]
    else:
        membraneblayer_list = []

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
    # 3D mesh algorithm (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=FrontalHex, 7=MMG3D, 9=R-tree), default = 1
    #raw_code(['Mesh.Algorithm3D = 1'])

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
