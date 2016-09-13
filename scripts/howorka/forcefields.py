"tools to save/load/interpolate/modify forcefields given by node-value data"

import numpy as np
import math
from itertools import product
import matplotlib.delaunay as dln
import dolfin
import nanopores
from nanopores.models import Howorka
from nanopores.tools import fields
import forcefield2D

# these params at least should be given to the force field computation
default = dict(
    Qmol = -1.,
    rMolecule = 0.5,
    Rx = 12.,
    Ry = 12.,
)
    
def forcefieldS1(overwrite=False, implicit=False, **params):
    # TODO: return right force field depending only on params (+ name)
    # as a first approximation, we will use:
    params1 = dict(default)
    params1.update(params)
    if implicit:           
        return forcefieldS1_implicit(overwrite, **params1)
    else:
        return forcefieldS1_explicit(overwrite, **params1)
        
def F_geo_phys(overwrite=False, implicit=False, **params):
    "returns force field and corresponding Geometry, Physics"
    params1 = dict(default)
    params1.update(params)
    if implicit:           
        F = forcefieldS1_implicit(overwrite, **params1)
    else:
        F = forcefieldS1_explicit(overwrite, **params1)
        
    mesh = F.function_space().mesh()
    geo, phys = Howorka.setup2D(mesh=mesh, z0=None, **params1)
    return F, geo, phys

# TODO: extend by zero on arbitrary size domain
def forcefieldS1_implicit(overwrite=False, **params):
    F, Fel, Fdrag, mesh, p = forcefield2D.maybe_calculate(overwrite, **params)
    return F

# TODO: "maybe_calculate" functionality could be implemented with decorator!
def forcefieldS1_explicit(overwrite=False, **params):
    # maybe interpolate and save force field
    if not overwrite and fields.exists(NAME, **params):
        print "Existing force field interpolation found."
    else:
        print "Interpolating 3D force field."
        interpolate_forcefieldS1_explicit(**params)
    
    # load force field
    FNAME = fields.get_entry(NAME, "FNAME", **params)
    forces, mesh, _ = nanopores.load_vector_functions(str(FNAME))
    F = forces["F"]
    return F
    
def load_forcefield_explicit(**params):
    FNAME = fields.get_entry(NAME, "FNAME", **params)
    forces, mesh, params = nanopores.load_vector_functions(str(FNAME))
    return forces["F"], forces["Fel"], forces["Fdrag"], mesh, params

NAME = "force2Dexp"
# TODO: put additional points also in point generation file
#       that is the only remaining ungeneric part of the routine!
def interpolate_forcefieldS1_explicit(**params):
    # --- load force field and points
    data = fields.get_fields("force3D", **params)
    X = data["x"]
    F = data["F"]
    Fexp = [[1e-12*f[0], 1e-12*f[2]] for f in F]
    Fimpl = forcefieldS1_implicit(**params)
    
    # --- load geometry params
    # TODO: not very general
    xparams = fields.load_file("xforce")["params"]
    Rx0, Ry0 = xparams["Rx"], xparams["Ry"]
    Ry0 = 7. # FIXME bad hack
    r = params["rMolecule"]
    r0, l0, r1, l1 = (Howorka.params_geo.r0, Howorka.params_geo.l0,
                     Howorka.params_geo.r1, Howorka.params_geo.l1)
    Rx, Ry = params["Rx"], params["Ry"]
        
    # --- piece together data TODO
    # additional points where implicit model will be used
    hadd = 0.8
    Xadd = uniform_grid([r1+r, Rx], [-Ry, -l1/2-r], hadd) + \
           uniform_grid([r1+r, Rx], [l1/2+r, Ry], hadd) + \
           uniform_grid([0, r1], [Ry0, Ry], hadd) + \
           uniform_grid([0, r1], [-Ry, -Ry0], hadd)
   # points where force is zero
    hdna = 0.4
    Xdna = uniform_grid([r0-r+.01, r1+r-.01], [-l0/2-r, l0/2+r], hdna) + \
           uniform_grid([r1+r, Rx], [-l1/2-r, l1/2+r], hdna)
    Fexpadd = [Fimpl(x) for x in Xadd]
    X += Xadd
    Fexp += Fexpadd
    Fimp = [Fimpl(x0) for x0 in X]
    X += Xdna
    Fdna = [[0.,0.] for x in Xdna]
    Fexp += Fdna
    Fimp += Fdna
    Fexp = np.array(Fexp)
    Fimp = np.array(Fimp)
    x = np.array([t[0] for t in X])
    y = np.array([t[1] for t in X])
    # overwrite explicit calculation far from pore # FIXME bad hack
    yfar = abs(y) > Ry0
    Fexp[yfar,:] = Fimp[yfar,:]
    
    # --- obtain S1 functions
    mesh = Fimpl.function_space().mesh()
    F = data_to_S1(x, y, mesh, Fexp=Fexp)
    
    # --- save functions
    uid = fields._unique_id()
    FNAME = NAME + uid
    nanopores.save_functions(FNAME, mesh, meta=params, F=F)
    fields.save_entries(NAME, params, FNAME=FNAME)
    fields.update()

def uniform_grid(xi, yi, h):
    Lx = xi[1]-xi[0]
    Ly = yi[1]-yi[0]
    m = int(math.ceil(Lx/h))
    n = int(math.ceil(Ly/h))
    x = np.linspace(xi[0], xi[1], m)
    y = np.linspace(yi[0], yi[1], n)
    X = list(product(x, y))
    return X

def data_to_S1(x, y, mesh, **values):
    "2D data clouds of vector valued functions => dolfin S1 functions"
    trid = dln.Triangulation(x, y)
    #mesh = nanopores.RectangleMesh([0, -Ry], [Rx, Ry], Nx, Ny)
    
    functions = []
    for F in values.values():    
        Fi = [None]*2
        for i in (0,1):
            intp = trid.nn_interpolator(F[:,i])
            intp2 = lambda x: intp([x[0]], [x[1]])
            Fi[i] = lambda_to_S1(intp2, mesh)
    
        V = dolfin.VectorFunctionSpace(mesh, "CG", 1)
        FF = dolfin.Function(V)
        dolfin.assign(FF, [Fi[0], Fi[1]])
        functions.append(FF)
    
    if len(functions) == 1:
        return functions[0]
    else:
        return tuple(functions)
    
def lambda_to_S1(f, mesh, dim=1):
    V = dolfin.FunctionSpace(mesh, "CG", 1)
    if dim>1:
        V = dolfin.MixedFunctionSpace([V]*dim)
    f1 = dolfin.Function(V)
    value_shape = () if dim==1 else (dim,)
    class expr(dolfin.Expression):
        def eval(self, value, x):
            value[:] = f(x)
        def value_shape(self):
            return value_shape
    f1.interpolate(expr())
    return f1
    
if __name__ == "__main__":
    from plot_forcefield import porestreamlines
    F = forcefieldS1(implicit=False, **default)
    Fimp = forcefieldS1(implicit=True, **default)
    porestreamlines(Howorka.polygon(), 6., 8., F=F, Fimp=Fimp)
    nanopores.showplots()
            