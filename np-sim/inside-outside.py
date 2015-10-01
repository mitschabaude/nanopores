from nanopores.geometries.H_geo.params_geo import *
from dolfin import *

def subdomain_list(**params):
    globals().update(params)
    return [AutoSubDomain(lambda _:True), MoleculeNotInFluid(), MoleculeInFluid()]
    
def boundaries_list(**params):
    globals().update(params)
    return []
    
synonymes = {
    "unclear":"autosubdomain",
}

r = rMolecule*1.1
    
def near_line(x0, X0): # line x[0] == X_0
    return between(x0, (X0-r,X0+r))
    
def near_point(x, X): # point X == [X[0],X[1]]
    d = sqrt((x[0]-X[0])**2 + (x[1]-X[1])**2)
    return between(d, (0,r))
    
def molecule_in_fluid(x):
    in_dna = between(x[0], (r0,r1)) and between(abs(x[1]), (0,0.5*l0))
    in_mem = between(x[0], (r1,Rx)) and between(abs(x[1]), (0,0.5*l1))
        
    near_dnab = ( (near_line(x[0], r0) or near_line(x[0], r1)) and between(abs(x[1]), (0,0.5*l0)) ) \
             or (  near_line(abs(x[1]), 0.5*l0)                and between(x[0], (r0,r1))         )
    near_memb = ( (near_line(x[0], r1) or near_line(x[0], Rx)) and between(abs(x[1]), (0,0.5*l1)) ) \
             or (  near_line(abs(x[1]), 0.5*l1)                and between(x[0], (r1,Rx))         )         
    near_boxb = near_line(x[0], Rx) or near_line(abs(x[1]), Ry)
    
    near_fluidb = near_dnab or near_memb #or near_boxb
                 
    corners = [[r0,0.5*l0], [r1,0.5*l0], [r1,-0.5*l0], [r0,-0.5*l0]]
    near_dna_corners = any([near_point(x, X) for X in corners])
    return not(in_dna or in_mem or near_fluidb or near_dna_corners)
    
    
class MoleculeInFluid(SubDomain):
    def inside(self, x, on_boundary):
        return molecule_in_fluid(x)
        
class MoleculeNotInFluid(SubDomain):
    def inside(self, x, on_boundary):
        return not molecule_in_fluid(x)
        

