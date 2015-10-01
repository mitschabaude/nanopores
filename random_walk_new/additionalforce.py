import math
from math import sqrt, pi
import numpy as np

# bad hack to allow arbitrary nm in params_geo
from nanopores import import_vars
params = import_vars("nanopores.W_3D_geo.params_geo")
for x in params:
    exec("%s = %s*1e0/%s" %(x, params[x], params['nm']))
def radius(x,y):
    return sqrt(x**2+y**2)

def additionalforce(x,y,z,delta=0,pitch=1):
    height = (lsin + lau + lsam) * nm
    r1 = height / math.tan(angle * pi / 180.) + r0 * nm
    k = math.tan(angle * pi / 180.)
    R = radius(x,y)
    if not(-height/2. <= z and z <= height/2.):
        return [0,0,0]
    if R < z + k*r0*nm + height/2. - (rMolecule*nm + delta)/(math.cos(pi/2.-angle*pi/180.)):
        return [0,0,0]
    coeff = (R - z - k*r0*nm - height/2. + (rMolecule*nm + delta)/(math.cos(pi/2.-angle*pi/180.)))*pitch/(sqrt((k**2+1)*(R**2)))
    return [-x*coeff,-y*coeff,R*coeff/k]



