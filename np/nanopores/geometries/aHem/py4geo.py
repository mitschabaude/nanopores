import numpy, math
from params_geo import *
import nanopores.py4gmsh.basic
import nanopores.py4gmsh.extra
from nanopores.py4gmsh import *
from warnings import warn

def get_geo(x0 = None, crosssections = True, **params):
    """
    writes a 3d geo file for an extruded axissymmetric geometry
    for Howorka 'Self-Assembled aHem that spans lipid bilayers'
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




    lt = -l2-(l0+l1-l2-lm)/2.
    lb = -l2-(l0+l1-l2-lm)/2.-lm


    X_Fluid_ext = numpy.array([[0.0, 0.0, l3],
                               [R, 0.0, l3],
                               [R, 0.0, -5.41],
                               [R, 0.0, -9.84],
                               [R, 0.0, -l0-l1-l4],
                               [0.0, 0.0, -l0-l1-l4]])
    X_Fluid_ctr = numpy.array([[0.0, 0.0, -10.21],#########################
                               [0.0, 0.0, -6.25],
                               [0.0, 0.0, -3.22],
                               [0.0, 0.0, 0.0]])
    X_aHem = numpy.array([[2.16, 0.0, 0.0],
                          [2.77, 0.0, -0.19],
                          [3.24, 0.0, -0.1 ],
                          [3.59, 0.0, -0.1 ],
                          [3.83, 0.0, -0.35],
                          [3.84, 0.0, -0.8 ],
                          [3.67, 0.0, -1.34],
                          [3.73, 0.0, -1.96],
                          [3.93, 0.0, -2.31],
                          [4.23, 0.0, -2.67],
                          [4.44, 0.0, -2.81],
                          [4.33, 0.0, -3.25],
                          [4.01, 0.0, -3.5 ],
                          [3.99, 0.0, -3.67],
                          [4.11, 0.0, -3.94],
                          [4.39, 0.0, -4.12],
                          [4.44, 0.0, -4.52],
                          [4.73, 0.0, -4.86],
                          [4.96, 0.0, -5.41],
                          [4.89, 0.0, -5.87],
                          [4.63, 0.0, -6.44],
                          [4.43, 0.0, -6.96],
                          [4.07, 0.0, -7.32],
                          [3.71, 0.0, -7.51],
                          [3.46, 0.0, -7.36],
                          [3.41, 0.0, -7.1 ],
                          [3.31, 0.0, -6.9 ],
                          [3.04, 0.0, -6.87],
                          [2.73, 0.0, -6.73],
                          [2.41, 0.0, -6.6 ],
                          [2.17, 0.0, -6.41],
                          [1.97, 0.0, -6.23],
                          [1.84, 0.0, -6.03],
                          [1.76, 0.0, -5.87],
                          [1.54, 0.0, -5.87],
                          [1.4 , 0.0, -5.96],
                          [1.31, 0.0, -6.16],
                          [1.39, 0.0, -6.57],
                          [1.6 , 0.0, -6.81],
                          [1.71, 0.0, -7.09],
                          [1.76, 0.0, -7.32],
                          [1.67, 0.0, -7.65],
                          [1.44, 0.0, -7.81],
                          [1.49, 0.0, -8.06],
                          [1.56, 0.0, -8.36],
                          [1.44, 0.0, -8.61],
                          [1.43, 0.0, -8.79],
                          [1.44, 0.0, -9.1 ],
                          [1.6 , 0.0, -9.48],
                          [1.74, 0.0, -9.84],
                          [1.63, 0.0, -10.0],
                          [1.47, 0.0, -10.19],
                          [1.26, 0.0, -10.21],
                          [1.07, 0.0, -10.05],
                          [1.03, 0.0, -9.76],
                          [1.09, 0.0, -9.44],
                          [1.07, 0.0, -9.02],
                          [0.86, 0.0, -8.79],
                          [0.64, 0.0, -8.68],
                          [0.63, 0.0, -8.36],
                          [0.8 , 0.0, -8.22],
                          [0.81, 0.0, -7.93],
                          [0.89, 0.0, -7.71],
                          [1.04, 0.0, -7.51],
                          [1.1 , 0.0, -7.25],
                          [0.91, 0.0, -7.02],
                          [0.91, 0.0, -6.76],
                          [0.91, 0.0, -6.48],
                          [0.69, 0.0, -6.25],
                          [0.69, 0.0, -6.  ],
                          [0.66, 0.0, -5.68],
                          [0.59, 0.0, -5.36],
                          [0.53, 0.0, -5.12],
                          [0.54, 0.0, -4.92],
                          [0.79, 0.0, -4.84],
                          [1.03, 0.0, -4.89],
                          [1.21, 0.0, -4.7 ],
                          [1.36, 0.0, -4.42],
                          [1.49, 0.0, -4.16],
                          [1.66, 0.0, -3.92],
                          [1.66, 0.0, -3.7 ],
                          [1.8 , 0.0, -3.41],
                          [2.  , 0.0, -3.22],
                          [1.91, 0.0, -2.93],
                          [1.8 , 0.0, -2.71],
                          [1.56, 0.0, -2.55],
                          [1.46, 0.0, -2.38],
                          [1.3 , 0.0, -2.19],
                          [1.21, 0.0, -1.93],
                          [1.09, 0.0, -1.64],
                          [0.9 , 0.0, -1.45],
                          [0.8 , 0.0, -1.28],
                          [0.84, 0.0, -1.  ],
                          [1.  , 0.0, -0.8 ],
                          [1.26, 0.0, -0.64],
                          [1.7 , 0.0, -0.31]])



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

    #Anchor Points on aHem for membran (index)
    ap1 = 18
    ap2 = 49
    apdiff=ap2-ap1

    e_Membrane = [Line(p_aHem[ap1],p_Fluid[2]), Line(p_Fluid[3], p_aHem[ap2])]


    edges_to_rot = [e_Fluid[0:5], e_aHem, e_Membrane]


    geo_cs_str = "no crosssectional surface"
    if crosssections:
        Comment(' integrate crosssectional lines in fluid and check if molecule intersects lines')
        e_CrossS = [Line(p_aHem[52], p_Fluid[6]), Line(p_aHem[68], p_Fluid[7]), Line(p_aHem[82], p_Fluid[8]), Line(p_aHem[0], p_Fluid[9])]
        cs_pop_i = None
        # check if molecule is near pore
        if x0 is not None and (x0[0]**2 + x0[1]**2 <= r0**2):
            # check z coordinate of molecule
            if abs(x0[2]) < rMolecule:
                geo_cs_str = "top crosssection"
                cs_pop_i = -1
            elif abs(x0[2] + 3.22) < rMolecule:
                geo_cs_str = "center top crosssection"
                cs_pop_i = 2
            elif abs(x0[2] + 3.25) < rMolecule:
                geo_cs_str = "center bottom crosssection"
                cs_pop_i = 1
            elif abs(x0[2] + 10.21) < rMolecule:
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
    
    tostr = lambda l: "{%s}"%(",".join(l),)
    ps_Fluid = PhysicalSurface(tostr(surfs_Fluid),'3') #Physical Surface Fluid
    
    surfs_Fluid_aHem = surfs[1][:]
    for index in range(apdiff):
        del surfs_Fluid_aHem[ap1::n_e_i[1]-index]  # deletes membrane
    [surfs_Fluid.append(s) for s in surfs_Fluid_aHem]
    [surfs_Fluid.append(s) for s in surfs[2]]
    sl_Fluid = SurfaceLoop(surfs_Fluid)

    ps_aHem = PhysicalSurface(tostr(surfs_Fluid_aHem),'1') #Physical Surface aHem

    sl_aHem = SurfaceLoop(surfs[1])
    vol_aHem = Volume(sl_aHem)

    surfs_Membrane = surfs[0][2::n_e_i[0]]
    for index in range(apdiff):
        [surfs_Membrane.append(s) for s in surfs[1][ap1+index::n_e_i[1]]]
    [surfs_Membrane.append(s) for s in surfs[2]]
    sl_Membrane = SurfaceLoop(surfs_Membrane)
    vol_Membrane = Volume(sl_Membrane)
    
    surfs_Membrane_ps = surfs[2]
    ps_Membrane = PhysicalSurface(tostr(surfs_Membrane_ps),'2') #Physical Surface Membrane


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
        
    PhysicalVolume(vol_Fluid, 1)
    PhysicalVolume(vol_Membrane, 2)
    PhysicalVolume(vol_aHem, 3)
    
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
                #"crosssections": crosssections,
                #"Number of crosssections": len(e_CrossS),
                #"Total number of crossections": 4,
                #"molecule crosses": geo_cs_str,
                #"popped crossection index": cs_pop_i,
                #"cs_pop_i": cs_pop_i,
                "Typical length scale on aHem": lcCenter,
                "geo_code": get_code(),
            }
    return geo_dict


# -----
if __name__ == '__main__':
    print(get_geo())
    print('\n - This is the sample code for the geo file')
