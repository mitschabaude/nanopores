'''
 provide default values and functions for calculating physical parameters
 for a specific physical set-up: nanopore with molecule inside
'''
import dolfin
from nanopores.physics.default import *
# TODO: we could easily "subclass" a more general physics module by importing * from it
#       e.g. default -> pore -> pore_molecule

# 1. -- default values for direct parameters

bV0 = None # voltage bias across pore [V] (None for no enforced bias)

Qmol = 0.*qq # charge on one molecule [C] = [A*s]
Membraneqs = -0.0 # membrane surface charge density [C/m**2]
DNAqsPure = -qq/nm**2 # = -0.16 # DNA surface charge density [C/m**2]
dnaqsdamp = 1. # DNA charge damping
SiNqs = -0.022
SAMqs = -0.078
ahemqs = 0.

T = 293 # temperature [K]
bulkcon = 300. # bulk concentration of ions [mol/m**3]
D = 1.9e-9  # diffusivity [m^2/s]

uppermbias = None # upper/lower concentration bias of negative/positive ions [mol/m**3]
lowermbias = None
upperpbias = None
lowerpbias = None

rpermPore = rpermw
rpermProtein = 2. # TODO ?????
rDPore = 0.5
stokesdampPore = 1.0

rTarget = 0.5*nm
qTarget = -qq
rDtargetPore = 1.
DtargetBulk = lambda: kB*T/(6*dolfin.pi*eta*rTarget) # Stokes Law
DtargetPore = lambda: kB*T/(6*dolfin.pi*eta*rTarget)*rDtargetPore
Dtarget = { # diffusivity of target molecule relative to bulk
    "bulkfluid": "DtargetBulk",
    "pore": "DtargetPore",
    "solid": 0.,
}

bulkconFluo = 0. # bulk concentration of fluorophore [mol/m**3]
hReservoir = 0.01 # height of cylindrical upper reservoir
applylowerqs = False
couplebVtoQmol = False
exactMqv = False
adaptMqv = True

# FIXME: add a list of functionals that can generically be used by PNPS or any system as its results
# --> because which functionals are calculated depends on the Physics of the problem!
# --> this could recognize which functionals make sense according to PDESystem.functions.keys()

# 2. -- derived parameters depending on other parameters and/or geometry
#    -- these are FUNCTIONS with keyword arguments corresponding to direct parameters

DNAqs = lambda: DNAqsPure*dnaqsdamp
permPore = lambda: eperm*rpermPore
permProtein = lambda: eperm*rpermProtein
kT = lambda: kB*T
UT = lambda: kB*T/qq
mu = lambda: D*qq/(kB*T) # mobility [m^2/Vs]
cFarad = lambda: qq*mol  # Faraday constant [C/mol]
debye = lambda: dolfin.sqrt(rpermw*eperm*kB*T/qq**2/2/mol/bulkcon)  # debye length [m]
bulkconduct = lambda: 2.*bulkcon*qq*mol*D*qq/(kB*T)  # 2*c0*cFarad*mu # electrolyte bulk conductivity [S/m]
lowerqs = lambda: (-Qmol*bulkconFluo*mol*hReservoir/2. if applylowerqs else 0.)

bV = lambda: -bV0*Qmol/qq if couplebVtoQmol else None

def Moleculeqs(geo, Qmol): # Molecule surface charge density [C/m**2]
    try:
        lscale = geo.parameter("nm")/nm
        scale = dolfin.Constant(1.0/lscale**2)
        r = dolfin.Expression("2*pi*x[0]")*scale if geo.params["dim"] == 2 else scale
        MolArea = dolfin.assemble(r('+')*geo.dS("moleculeb"))
        return Qmol/MolArea if MolArea > 0. else 0.
    except Exception:
        return None

def Moleculeqv(geo, Qmol, exactMqv, adaptMqv, lscale): # Molecule volume charge density [C/m**3]
    if exactMqv:
        r = geo.parameter("rMolecule")/lscale
        vol = 4.*dolfin.pi/3.*r**3
        print "Molqv: ", Qmol/vol
        return Qmol/vol if vol > 0. else 0.
        #return -327926363.681 #-305959545.378
    elif adaptMqv:
        scale = dolfin.Constant(1.0/lscale**3)
        r = dolfin.Expression("2*pi*x[0]")*scale if geo.params["dim"] == 2 else scale
        def compute(geo):
            vol = dolfin.assemble(r*geo.dx("molecule"))
            return Qmol/vol if vol > 0. else 0.
        const = geo.constant("Moleculeqv", compute)
        print "Molqv: ", geo.constants["Moleculeqv"].value
        return const
    try:
        lscale = geo.parameter("nm")/nm
        scale = dolfin.Constant(1.0/lscale**3)
        r = dolfin.Expression("2*pi*x[0]")*scale if geo.params["dim"] == 2 else scale
        MolVol = dolfin.assemble(r*geo.dx("molecule"))
        return Qmol/MolVol if MolVol > 0. else 0.
    except Exception:
        return None
        
def DNAArea(geo, lscale):
    # this is for Howorka DNA nanopores
    h = geo.params["l0"] # height of DNA
    h1 = geo.params["l1"] # height of membrane
    r0 = geo.params["r0"] # inner radius
    r1 = geo.params["r1"] # outer radius
    pi = dolfin.pi
    geo.constant("DNAArea0", lambda geo: 2.*pi*(h*r0 + (h-h1)*r1)/lscale**2)
    return 2.*pi*(h*r0 + (h-h1)*r1)/lscale**2
    
def QDNA(DNAqs, DNAArea):
    # total charge on DNA if boundary is exactly cylindrical, i.e. in 2D
    return DNAqs*DNAArea
        
def DNAqsHoworka(geo, DNAqs, QDNA, dim, r2pi, lscale):
    if dim == 2:
        return DNAqs
    # calculate total charge on DNA if boundary is exactly cylindrical
    scale = dolfin.Constant(1.0/lscale**2)
    def compute(geo):
        area = dolfin.assemble(dolfin.avg(scale*r2pi)*geo.dS("chargeddnab"))
        return QDNA/area if area > 0. else 0.
    geo.constant("DNAArea", lambda geo: dolfin.assemble(dolfin.avg(scale*r2pi)*geo.dS("chargeddnab")))
    
    return geo.constant("DNAqs", compute)

# 3. -- piecewise maps connecting physical domains to parameters
#    -- these are DICTIONARIES of the form {"domain": "parameter_name"} OR
#                                          {"domain": parameter}
permittivity.update(
    bulkfluid = eperm*rpermw,
    pore = "permPore",
    protein = "permProtein",
)

# determines how Molecule charge is implemented
smearMolCharge = True # note: this is no parameter

surfcharge = { # surface charge densities for Neumann RHS
    "moleculeb": (0. if smearMolCharge else "Moleculeqs"),
    "chargedmembraneb": "Membraneqs",
    "chargeddnab": "DNAqsHoworka", #"DNAqs",
    "chargedsinb": "SiNqs",
    "chargedsamb": "SAMqs",
    "lowerb": "lowerqs",
    "ahemb": "ahemqs",
}

volcharge = dict( # volume charges for RHS
    default = 0.,
    molecule = ("Moleculeqv" if smearMolCharge else 0.),
)
    
charge = {"volcharge":volcharge, "surfcharge":surfcharge}

diffusion_factor = { # diffusivity of ions relative to bulk
    "bulkfluid": 1.,
    "pore": "rDPore",
    "solid": 0.,
}

stokes_damp = { # damps stokes terms in current densities
    "bulkfluid": 1.,
    "pore": "stokesdampPore",
    "solid": 0.,
}

initial_ions = {
    "fluid": "bulkcon",
    "solid": "bulkcon",
}

def r2pi(geo, dim):
    return dolfin.Expression("2*pi*x[0]") if dim == 2 else dolfin.Constant(1.0)


# TODO: i don't know yet how we should handle functionals
'''
Scaling for functionals:
in our package all volume forms are in fact scaled with
1./lscale**3, thus this is the reference for all forms,
this means surfaces and grad need to be scaled by lscale
'''
def Fbare(geo, r2pi, Moleculeqs, Moleculeqv, grad, lscale):
    try: # to make geo not necessary
        if len(geo.physicalboundary("moleculeb"))==0:
            return
        def Fbaresurf(v, i):
            dS = geo.dS("moleculeb")
            n = dolfin.FacetNormal(geo.mesh)
            #def tang(x):
            #    return x - dolfin.inner(x,n)*n
            #def gradT(v): # tangential gradient
            #    return tang(dolfin.grad(v))
            #return dolfin.Constant(Moleculeqs)*(-r*gradT(v)[i])('-')*dS
            return dolfin.Constant(Moleculeqs) * (-r2pi*grad(v)[i])('-')*dS
        def Fbarevol(v, i):
            dx = geo.dx("molecule")
            scale = dolfin.Constant(lscale**(-3))
            return Moleculeqv * (-r2pi*grad(v)[i])*dx

        if geo.params["x0"]:
            Fbare0 = Fbarevol if smearMolCharge else Fbaresurf
        else: Fbare0 = None

        return Fbare0
    except:
        return None

def FbareE(geo, r2pi):
    try:
        # for gradient recovery, or testing with manifactured E field
        def Fbaresurf(E, i):
            dS = geo.dS("moleculeb")
            return dolfin.Constant(Moleculeqs) * (-r2pi*E[i])('-') * lscale*dS
        def Fbarevol(E, i):
            dx = geo.dx("molecule")
            return dolfin.Constant(Moleculeqv) * (-r2pi*E[i]) * dx

        if geo.params["x0"]:
            Fbare0 = Fbarevol if smearMolCharge else Fbaresurf
        else: Fbare0 = None

        return Fbare0
    except:
        return None

def CurrentPB(geo, r2pi, bulkcon, mu, rDPore, bV, UT, lscale, cFarad):
    try:
        # current approximated by linear functional of potential (not really)
        bV0 = 0.1 #bV if bV is not None else 0.1
        def J0(v):
            L = geo.params["l0"]/lscale
            E = bV0/L
            Jz = 2*cFarad*bulkcon*mu*rDPore*v/UT*E* r2pi/L*geo.dx("pore")
            return Jz
        return J0
    except:
        return None

def CurrentPBdrift(geo, r2pi, bulkcon, mu, rDPore, UT, lscale, cFarad):
    try:
        def Jzdrift(bV):
        #def Jzdrift(v):
            # drift current approximated by PB
            L = geo.params["l0"]/lscale
            dim = geo.params["dim"]
            #return dolfin.assemble(2*cFarad*bulkcon*mu*rDPore*(-grad(v)) * r2pi/L*geo.dx("pore"))
            return dolfin.assemble(2*cFarad*bulkcon*mu*rDPore*bV/L * r2pi/L*geo.dx("pore"))
        return Jzdrift
    except:
        return None
        
def Feff(geo, grad, qTarget, rTarget):
    def Feff0(v, u):
        E = -grad(v)
        pi = 3.141592653589793
        Fel = dolfin.Constant(qTarget)*E
        Fdrag = dolfin.Constant(6*pi*eta*rTarget)*u
        F = Fel + Fdrag
        return F
    return Feff0
