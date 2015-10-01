import numpy
from params_geo import *
import nanopores.py4gmsh.basic
import nanopores.py4gmsh.extra
from nanopores.py4gmsh import *

def get_geo(**params):
    """
    writes a 2d geo file for an axisymmetric geometry for particle in cylinder
    """
    basic = reload(nanopores.py4gmsh.basic)
    extra = reload(nanopores.py4gmsh.extra)
    globals().update(params)
    
    lcMolecule = r/10
    lcFluid = min(min(R,l/2)/5, max(R,l/2)/60)

    X_Fluid = numpy.array([[0, l/2, 0],
                           [R, l/2, 0],
                           [R,-l/2, 0],
                           [0,-l/2, 0]])
                           
    X_Molecule = numpy.array([[0, z0 - r, 0],
                              [0, z0    , 0],
                              [0, z0 + r, 0]])

    p_Fluid = [Point(x, lcFluid) for x in X_Fluid]
    p_Molecule = [Point(x, lcMolecule) for x in X_Molecule]
    
    p_Fluid.append(p_Molecule[0])
    p_Fluid.append(p_Molecule[2])

    c_Molecule = Circle(p_Molecule)
    e_Molecule = Line(p_Molecule[2], p_Molecule[0])
    ll_Molecule = LineLoop([c_Molecule, e_Molecule])
    s_Molecule = PlaneSurface(ll_Molecule)

    e_Fluid = [Line(p_Fluid[k], p_Fluid[k+1]) for k in range(len(p_Fluid)-1)]
    e_Fluid.append(Line(p_Fluid[-1], p_Fluid[0]))
    e_Fluid.pop(len(e_Fluid)-2)
    e_Fluid.insert(len(e_Fluid)-1, c_Molecule)
    ll_Fluid = LineLoop(e_Fluid)
    s_Fluid = PlaneSurface(ll_Fluid)
    
    box1 = BoxField(lcFluid, lcFluid, 0.0, R, -l/2, l/2)
    box2 = BoxField(lcMolecule, lcFluid, -1e-15, 6*r, z0-6*r, z0+6*r)
    field_list = [box1,box2]
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
                "Typical length scale": lcFluid,
                "geo_code": get_code(),}
    return geo_dict


# -----
if __name__ == '__main__':
    print(get_geo())
    print('\n - This is the sample code for the geo file')
