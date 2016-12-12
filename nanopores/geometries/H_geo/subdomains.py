"""
mark domains/boundaries with dolfin MeshFunctions
"""

"""
ueberlegungen bzgl. abstraktion:

subdomain numbers werden gebraucht fuer:
-) bilinearformen (masse) --> als tupel
   -> jede variable in (v,ionen,wasser) ist mit einem tupel \subset (0,1,2,..) assoziiert
   -> diese assoziation ist eine geometrische eigenschaft, gehoert daher in geo file
-) fuer potential: permittivity
-) fuer ionen: versch. diffusionskonstanten
-) spaeter fuer wasser vielleicht etwas aehnliches, bei unified continuum mechanics (??)
-) charges, d.h. neumann boundaries fuer potential
   -> moeglicherweise sowohl auf interfaces als auch boundaries
   -> bisher aber nur auf interfaces --> also mass dS
   -> tupel \subset (0,1,2,..) im geo file passend zu boundaries_list()
      UND tupel mit den entsprechenden charges
   -> sowas wie

	Lv = 0.0*dS
	for c in chargedinterfaces
		Lv = Lv + Constant(charge(c))*dS(c)

-) dirichlet boundaries + randwerte

es waere schoen das alles in einer klasse zu haben und nicht alles erst "getten" zu muessen
"""

from dolfin import *
from math import sqrt, pow
from . import params_geo

synonymes = {
    "bulkfluid": {"bulkfluidtop", "bulkfluidbottom"},
    "fluid":{"bulkfluid","pore"},
    "pore":{"poretop", "porecenter", "porebottom"},
    "lipid":"membrane",
    "solid":{"dna", "membrane", "molecule"},
    "ions":"fluid",
    "bulk":{"upperb","lowerb"}, #,"rightfluidb"},
    "chargeddnab":{"chargeddnainb","chargeddnaoutb"},
    "dnab":{"chargeddnab","unchargeddnab"},
    "noslip":{"dnab","membraneb"}, #"moleculeb"},
    "nopressure":{"upperb","lowerb"},
    #"charged":{"chargeddnab","moleculeb","membraneb"},
    "ground":"upperb",
    "bV":"lowerb",
    "chargedmembraneb":"membraneb",
}

def norm2(x, y):
    return sqrt(pow(x,2) + pow(y,2))

# lists containing subdomain classes, ordering is important: fluid first, molecule last
def subdomain_list(**params):
    tmp = vars(params_geo)
    params.update({key:tmp[key] for key in tmp  \
                   if (not key in params and not key.startswith("__"))})
    tolc = params["tolc"]
    try:
        x0 = params["x0"]
    except (KeyError, NameError):
        x0 = None

    # subdomains
    class BulkFluidTop(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] > -tolc  # other domains will overwrite

    class BulkFluidBottom(SubDomain):
        def inside(self, x, on_boundary):
            return x[1] < tolc  # other domains will overwrite

    class Molecule(SubDomain):
        def inside(self, x, on_boundary):
            if x0 is not None:
                return norm2(x[0], x[1]-x0[2]) <= params["rMolecule"] +tolc
            else:
                return False

    #class MoleculeHull(SubDomain):
    #    def inside(self, x, on_boundary):
    #        if x0 is not None:
    #            return norm2(x[0], x[1]-x0[2]) <= params["rMolecule"] + params["rMH"] +tolc
    #        else:
    #            return False

    class DNA(SubDomain):
        def inside(self, x, on_boundary):
            return (x[0] >= (params["r0"] - tolc) and x[0] <= (params["r1"] +tolc)  \
                    and (abs(x[1])-tolc) <= 0.5*params["l0"])

    class Membrane(SubDomain):
        def inside(self, x, on_boundary):
            return (x[0] >= (params["r1"] -tolc) and x[0] <= params["Rx"]  \
                    and abs(x[1]) -tolc <= params["l1"]/2 )

    # partion pore into three subdomains
    class PoreTop(SubDomain):
        def inside(self, x, on_boundary):
            return (between(x[1],(params["l1"]/2, params["l0"]/2)) and between(x[0], (0,params["r0"])))

    class PoreCenter(SubDomain):
        def inside(self, x, on_boundary):
            return (between(x[1],(-params["l1"]/2, params["l1"]/2)) and between(x[0], (0,params["r0"])))

    class PoreBottom(SubDomain):
        def inside(self, x, on_boundary):
            return (between(x[1],(-params["l0"]/2,-params["l1"]/2)) and between(x[0], (0,params["r0"])))

    return [BulkFluidTop(), BulkFluidBottom(), DNA(), Membrane(),
            PoreTop(), PoreCenter(), PoreBottom(), Molecule(),]

def boundaries_list(**params):
    tmp = vars(params_geo)
    params.update({key:tmp[key] for key in tmp  \
                   if (not key in params and not key.startswith("__"))})
    tolc = params["tolc"]
    try:
        x0 = params["x0"]
    except (KeyError, NameError):
        x0 = None

    # exterior fluid boundaries
    class UpperB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], params["Ry"])

    class LowerB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], -params["Ry"])

    class LeftFluidB(SubDomain):
        def inside(self, x, on_boundary):
            if x0 is not None:
                return on_boundary and near(x[0], 0) and  \
                    not between(x[1], (x0[2] - params["rMolecule"], x0[2] + params["rMolecule"]))
            else:
                return on_boundary and near(x[0], 0)

    class RightFluidB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], params["Rx"]) and abs(x[1]) >= params["l1"]/2 -tolc

    # DNA boundaries
    class ChargedDNAinB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (params["r0"] -tolc, params["r0"] +tolc))  \
                and between(abs(x[1]), ( -tolc, params["l0"]/2 +tolc))

    class ChargedDNAoutB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (params["r1"] -tolc, params["r1"] +tolc)) \
                and between(abs(x[1]), (params["l1"]/2 -tolc, params["l0"]/2 +tolc))

    class UnchargedDNAB(SubDomain):
        def inside(self, x, on_boundary):
            return ( between(x[0], (params["r0"] -tolc, params["r1"] +tolc)) \
                 and near(abs(x[1]), params["l0"]/2) )
            #return ( ( near(x[0], params["r0"]) and between( x[1], (-params["l1"]/2 -tolc, params["l1"]/2 + tolc) ) )  \
            #        or  ( between(x[0], (params["r0"] -tolc, params["r1"] +tolc)) and near(abs(x[1]), params["l0"]/2) ) )

    # Molecule boundaries
    class MoleculeB(SubDomain):
        def inside(self, x, on_boundary):
            if x0 is not None:
                return between(norm2(x[0], x[1]-x0[2]), (params["rMolecule"] -tolc, params["rMolecule"] +tolc) )
            else:
                return False

    # Membrane boundaries
    class MembraneB(SubDomain):
        def inside(self, x, on_boundary):
             return between(x[0], (params["r1"] -tolc, params["Rx"] +tolc)) and near(abs(x[1]), params["l1"]/2)

    # cross-section interfaces
    class CrossTop2D(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[1],params["l0"]/2) and between(x[0],(0,params["r0"])))

    class CrossCenterTop2D(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[1],params["l1"]/2) and between(x[0],(0,params["r0"])))

    class CrossCenterBottom2D(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[1],-params["l1"]/2) and between(x[0],(0,params["r0"])))

    class CrossBottom2D(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[1],-params["l0"]/2) and between(x[0],(0,params["r0"])))

    return [UpperB(), LowerB(), LeftFluidB(), RightFluidB(),
            ChargedDNAinB(), ChargedDNAoutB(), UnchargedDNAB(), MembraneB(),
            CrossTop2D(), CrossCenterTop2D(), CrossCenterBottom2D(), CrossBottom2D(), MoleculeB(),]


#the following code seems to be only needed for backward compatibility

# get MeshFunctions
# def get_subdomain(mesh):
#     subdomain = CellFunction("size_t", mesh, 0)
#     for i,sub in enumerate(subdomain_list()):
#         sub.mark(subdomain, i+1)
#     return subdomain

# def get_boundaries(mesh):
#     boundaries = FacetFunction("size_t", mesh, 0)
#     UpperB().mark(boundaries, 11)
#     LowerB().mark(boundaries, 12)
#     LeftFluidB().mark(boundaries, 13)
#     RightFluidB().mark(boundaries, 14)
#     ChargedDNAinB().mark(boundaries, 21)
#     ChargedDNAoutB().mark(boundaries, 22)
#     UnchargedDNAB().mark(boundaries, 23)
#     MembraneB().mark(boundaries, 31)
#     #MoleculeB(x0).mark(boundaries, 41, False)
#     return boundaries

# def get_porepartitions(mesh, x0):
#     crosssections = get_subdomain(mesh, x0)
#     PoreTop().mark(crosssections, 51, False)
#     PoreCenter().mark(crosssections, 52, False)
#     PoreBottom().mark(crosssections, 53, False)
#     Molecule(x0).mark(crosssections, 4)
#     return crosssections

# # get piecewise constant (pwc) function on subdomains (e.g. permittivity)
# # specified by list/array either as CellFunction or as DG0 Function
# class CellFunction2Expression(Expression):
#     def __init__(self, cellfun):
#         self.cellfun = cellfun
#     def eval_cell(self, values, x, cell):
#         values[0] = self.cellfun[cell.index]

# def get_pwc_cellfun(mesh,somearray):
#     cellfun = CellFunction("double", mesh)
#     for i,sub in enumerate(subdomain_list()):
#         sub.mark(cellfun, somearray[i])
#     return cellfun

# def get_pwc_DG(mesh,somearray):
#     cellfun = get_pwc_cellfun(mesh,somearray)
#     expr = CellFunction2Expression(cellfun)
#     dgfun = Function(FunctionSpace(mesh,"DG",0))
#     dgfun.interpolate(expr)
#     return dgfun

# def get_permittivity_DG(mesh):
#     return get_pwc_DG(mesh,perm)

# def get_diffusion_DG(mesh):
#     return get_pwc_DG(mesh,diff)
