"""
define dolfin SubDomains, dynamically depending on geometric imput parameters
"""

from dolfin import *
from .params_geo import *

# lists containing subdomain classes, ordering is important
def subdomain_list(**params):
    globals().update(params)
    #return [Oxide(), Core(), Fluid()]
    return [Oxide(), Core(), Fluid(), Dopants(xi)] if len(xi)>0 else \
           [Oxide(), Core(), Fluid()]
               
def boundaries_list(**params):
    globals().update(params)
    return [Upper(), Lower(), CoreB(), OxideB()]
    
#def mark_directly():
#    return [Dopants]   
    
synonymes = {
    "ions":"fluid",
    "solid":{"core", "oxide"},
    "sil":"core",
    "bulk":"upper",
    "bV":"lower",
    "ground":"upper",
    #"charged":"coreb",
    #"ground":"oxideb",
    #"biasedvoltage":{"upper","lower","coreb"},
    "charged":set(),
    "noslip":"oxideb",
    "nopressure":"upper",
    "impurity":"dopants",
    "wire":{"dopants", "core"}
}

# helper functions
# don't use these for negations

def in_rect_xz(x, C, w, h):
    # test whether x is inside x-z plane rectangle defined by
    # lower left corner C = (C_x, C_z), height h and width w
    return between(x[0], (C[0], C[0] + w)) and between(x[2], (C[1], C[1] + h))
    
def on_rect_xz(x, C, w, h):
    # test whether x is on x-z plane rectangle boundary defined by
    # lower left corner C = (C_x, C_z), height h and width w
    return ( between(x[2], (C[1], C[1] + h)) and (near(x[0], C[0]) or near(x[0], C[0] + w)) ) or \
           ( between(x[0], (C[0], C[0] + w)) and (near(x[2], C[1]) or near(x[2], C[1] + h)) )
           
def dist(x, x0):
    return sqrt(sum((t-t0)**2 for (t,t0) in zip(x,x0)))

def inside_dopant(x, x0):
    # test whether x is inside dopant with center x0 and radius r
    #return dist(x[::2], x0[::2]) <= r_eff + tolc
    return dist(x, x0) <= r_eff + tolc
    
def mark_cells_by_points(cellfunction, points, value):
    from numpy import array
    btree = cellfunction.mesh().bounding_box_tree()
    for x in points:
        cells = btree.compute_entity_collisions(Point(array(x)))
        if cells.size == 0:
            warning("The point %s does not lie inside mesh and was therefore not considered for marking subdomain." 
            %repr(x))
        else:
            cellfunction.set_value(cells[0], value)
    
    
# subdomains that are marked directly (without SubDomain)
# these must expect the following input arguments:
#    -) subdomains ... cellfunction to be marked
#    -) i ... the value used for marking
#    -) any other parameters as keyword arguments
def Dopants2(subdomains, i, xi=[], **params):
    mark_cells_by_points(subdomains, xi, i)
    
        
# subdomains
class Oxide(SubDomain):
    def inside(self, x, on_boundary):
        return True  # other domains have to overwrite

class Core(SubDomain):
    def inside(self, x, on_boundary):
        return in_rect_xz(x, (x_core, z_core), w_core, h_core)

class Fluid(SubDomain):
    def inside(self, x, on_boundary):
        return (x[2] >= h_side + h_wire - tolc) or \
               in_rect_xz(x, (0.0, h_side), w_side, h_wire) or \
               in_rect_xz(x, (w_side + w_wire, h_side), w_side, h_wire)
               
class Dopants(SubDomain):
    def __init__(self, xi):
        SubDomain.__init__(self)
        self.xi = xi
    def inside(self, x, on_boundary):
        return any(inside_dopant(x, x0) for x0 in self.xi)
        
               
# boundaries
class Upper(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] >= h - tolc

class Lower(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] <= tolc

class CoreB(SubDomain):
    def inside(self, x, on_boundary):
        return on_rect_xz(x, (x_core, z_core), w_core, h_core)
        
class OxideB(SubDomain): # fluid-oxide interface
    def inside(self, x, on_boundary):
        return (near(x[2], h_side) or on_rect_xz(x, (w_side, h_side), w_wire, h_wire)) and \
               not (near(x[2], h_side) and (w_side + tolc <= x[0] <= w_side + w_wire - tolc))

