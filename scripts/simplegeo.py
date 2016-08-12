""" example for how to use Geometry class with module """

from dolfin import *

mesh = UnitSquareMesh(1000,1000)
subdomains = CellFunction("size_t",mesh,0)
boundaries = FacetFunction("size_t",mesh,0)

DomainBoundary().mark(boundaries,1)
CompiledSubDomain("on_boundary && near(x[0],0)").mark(boundaries,2)
CompiledSubDomain("on_boundary && near(x[0],1)").mark(boundaries,3)

CompiledSubDomain("x[0] <= 0.3 || x[0] >= 0.8").mark(subdomains,1)

physical_boundary = {
    'left':(2,),
    'right':(3,),
    'other':(1,)
}

physical_domain = {
    'fluid':(1,),
    'dna':(0,)
}

synonymes = {
    'not_left':{'right','other'},
    'noslip':'other',
    'nopressure':'right',
    'inflow':'left',
}




