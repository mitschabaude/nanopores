import dolfin
from nanopores.physics.params_physical import *

def lscale(geo):
    # TODO: "lscale" is confusing since it is actually 1 over the length scale
    try:
        return geo.parameter("lscale")
    except:
        try:
            return geo.parameter("nm")/nm
        except:
            return 1e9

def grad(lscale):
    def grad0(u):
        return dolfin.Constant(lscale)*dolfin.nabla_grad(u)
    return grad0

def div(lscale):
    def div0(u):
        return dolfin.Constant(lscale)*dolfin.transpose(dolfin.nabla_div(u))
    return div0

def dim(geo):
    return geo.mesh.topology().dim()

cyl = False
def r2pi(cyl):
    return dolfin.Expression("2*pi*x[0]", degree=1) if cyl else dolfin.Constant(1.)

def invscale(lscale):
    return lambda i: dolfin.Constant(1./lscale**i)