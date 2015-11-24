import dolfin
from nanopores.physics.params_physical import *

def lscale(geo):
    try:
        return geo.parameter("nm")/nm
    except:
        return 1.0
def grad():
    def grad0(u):
        return lscale*dolfin.nabla_grad(u)
    return grad0
def div():
    def div0(u):
        return lscale*dolfin.transpose(dolfin.nabla_div(u))
    return div0

