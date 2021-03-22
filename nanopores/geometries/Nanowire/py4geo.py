import numpy
from .params_geo import *
import nanopores.py4gmsh.basic
import nanopores.py4gmsh.extra
from nanopores.py4gmsh import *
import importlib

def polygonchain(*sequence):
    # creates polygon chain from sequence of direction vectors
    # returns array of polygon points
    a = []
    x = numpy.array([0.0, 0.0, 0.0])
    for v in sequence:
        x = x + numpy.array(v)
        a.append(x)
    return numpy.array(a)
    
def zigzag(x0, *sequence, **kw):
    # creates polygon chain in x-z plane with alternating z/x steps
    # x0 .. starting point
    # i ... dimension index of first step (alternates between 0 and 2)
    a = [x0]
    i = kw.get("i",2)
    for s in sequence:
        v = [0,0,0]
        v[i] = s
        #print v
        a.append(v)
        i = 2 if i==0 else 0
    #print a
    return polygonchain(*a)
    
#def get_geo(xi = [], **params):
def get_geo(**params):
    """
    writes geo file for a 3d nanowire geometry
     
     /__________ /
    |           |/
    |   fluid   |/
    |   _____   |/
    |__|/sil/|__|/
    |///////////|/
    |// oxide //|/
    |///////////|/
    
    #INPUT: xi should be a list of 3-vectors of coordinates of
    #centers of balls to be placed inside silicon core.
    """
    importlib.reload(nanopores.py4gmsh.basic)
    importlib.reload(nanopores.py4gmsh.extra)
    
    # let **params overwrite params_geo
    globals().update(params)
    lcwire = lc * rellcwire
    
    # create points in x-z plane
    X_core = zigzag([x_core, 0, z_core], h_core, w_core, -h_core)
    X_oxide = zigzag([0, 0, h_side], w_side, h_wire, w_wire, -h_wire, w_side, -h_side, -w, i=0)
    X_top = zigzag([0, 0, h], w, i=0)
    
    p_core = [Point(x, lcwire) for x in X_core]
    p_oxide = [Point(x, lc) for x in X_oxide]
    p_top = ([Point(x, lc) for x in X_top])

    # create line loops and plane surfaces
    e_core = [Line(p_core[k], p_core[k+1]) for k in range(-1,3)]
    e_ctr = [Line(p_oxide[k], p_oxide[k+1]) for k in range(0,5)] # first 5 lines of p_oxide
    e_bottom = [Line(p_oxide[k], p_oxide[k+1]) for k in range(-3,0)] # last 3 lines
    e_top = [Line(p_oxide[0], p_top[0]),
             Line(p_top[0], p_top[1]),
             Line(p_top[1], p_oxide[5])]
    
    e_ctr_reverse = ["-%s"%s for s in e_ctr[::-1]]

    ll_core = LineLoop(e_core)
    ll_oxide = LineLoop(e_ctr + e_bottom)
    ll_fluid = LineLoop(e_top + e_ctr_reverse)

    s_core = PlaneSurface(ll_core)
    s_fluid = PlaneSurface(ll_fluid)
    # Workaround PlaneSurface accepts only one LineLoop
    s_oxide = PlaneSurface(",".join([ll_core, ll_oxide]))
    
    # create 3D geometry by translating all surfaces
    #raw_code(["Extrude {0, %s, 0} {\n  Surface{%s, %s, %s};\n}"
    #           % (l, s_core, s_fluid, s_oxide)])
    Extrude('Surface{%s, %s}' %(s_fluid, s_oxide),
             translation_axis = [0., l, 0.])
    ex = Extrude('Surface{%s}' %(s_core),
             translation_axis = [0., l, 0.])
    # ex[1] = box, ex[0] = back, ex[2,3,4,5] = sides
    #core_surfs = [s_core, ex[0]] + [ex[i] for i in range(2,6)]
    
    # to disable question dialogs
    raw_code(['General.ExpertMode = 1;'])

    # Meshing Algorithm: 2= ?, 5 = frontal (netgen)
    #raw_code(['Mesh.Algorithm3D = 5;'])

    geo_dict = {"gmsh mesh generating sript": __name__,
                "geo_code": get_code(),}
    return geo_dict


# -----
if __name__ == '__main__':
    print((get_geo()))
    print('\n - This is the sample code for the geo file')
