'''
 provide default values and functions for calculating physical parameters
 for a specific physical set-up: nanopore with DNA inside,
in this setup moleucle and DNA is the same!!!
'''

import dolfin
from nanopores.physics.params_physical import *
from warnings import warn

# 1. -- default values for direct parameters

bV = 0.1 # voltage bias across pore [V] (None for no enforced bias)

Membraneqs = -0.03 # membrane surface charge density [C/m**2]
DNAqsPure = -qq/nm**2 # = -0.16 # DNA surface charge density [C/m**2]
dnaqsdamp = 1. # DNA charge damping
SiNqs = -0.022
SAMqs = -0.078
Moleculeqv = 0.  # no molecule/Dna volume charges!

T = 293 # temperature [K]
bulkcon = 1e3 # bulk concentration of ions [mol/m**3]
D = 1.9e-9  # diffusivity [m^2/s]

uppermbias = None # upper/lower concentration bias of negative/positive ions [mol/m**3]
lowermbias = None
upperpbias = None
lowerpbias = None

rpermPore = rpermw
rDPore = 0.5
stokesdampPore = 1.0

# 2. -- derived parameters depending on other parameters and/or geometry
#    -- these are FUNCTIONS with keyword arguments corresponding to direct parameters

DNAqs = lambda: DNAqsPure*dnaqsdamp
permPore = lambda: eperm*rpermPore
kT = lambda: kB*T
UT = lambda: kB*T/qq
mu = lambda: D*qq/(kB*T) # mobility [m^2/Vs]
cFarad = lambda: qq*mol  # Faraday constant [C/mol]
debye = lambda: dolfin.sqrt(rpermw*eperm*kB*T/qq**2/2/mol/bulkcon)  # debye length [m]
bulkconduct = lambda: 2.*bulkcon*qq*mol*D*qq/(kB*T)  # electrolyte bulk conductivity [S/m]

Moleculeqs = DNAqs
Moleculeqv = 0.

# scaling hack:
# for scaled gradient, use phys.grad
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


# 3. -- piecewise maps connecting physical domains to parameters
#    -- these are DICTIONARIES of the form {"domain": "parameter_name"} OR
#                                          {"domain": parameter}

permittivity = {
    "bulkfluid": eperm*rpermw,
    "pore": "permPore",
    "dna": eperm*rpermDNA,
    "molecule": eperm*rpermDNA,
    "sin": eperm*rpermSiN,
    "lipid": eperm*rpermLipid,
    "au": eperm*rpermAu,
    "sam": eperm*rpermSAM,
}

smearMolCharge = False # note: this is no parameter

surfcharge = { # surface charge densities for Neumann RHS
    "chargedmembraneb": "Membraneqs",
    "chargeddnab": "DNAqs",
    "chargedsinb": "SiNqs",
    "chargedsamb": "SAMqs",
}

volcharge = {# volume charges for RHS
    "ions": 0.,
    "dna": 0.,
    "membrane": 0.,
    }

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
    "solid": 0.,
}

def r(geo):
    try:
        return dolfin.Expression("2*pi*x[0]") if geo.params["dim"] == 2 else dolfin.Constant(1.0)
    except:
        return None

# TODO: i don't know yet how we should handle functionals
'''
Scaling for functionals:
in our package all volume forms are in fact scaled with
1./lscale**3, thus this is the reference for all forms,
this means surfaces and grad need to be scaled by lscale
'''
def Fbare(geo):
    try: # to make geo not necessary
        def Fbaresurf(v, i):
            dS = geo.dS("chargeddnab")
            n = dolfin.FacetNormal(geo.mesh)
            #def tang(x):
            #    return x - dolfin.inner(x,n)*n
            #def gradT(v): # tangential gradient
            #    return tang(dolfin.grad(v))
            #return dolfin.Constant(Moleculeqs)*(-r*gradT(v)[i])('-')*dS
            return dolfin.Constant(DNAqs) * (-r*lscale*dolfin.grad(v)[i])('-') * lscale*dS
        def Fbarevol(v, i):
            dx = geo.dx("molecule")
            return dolfin.Constant(Moleculeqv) * (-r*lscale*dolfin.grad(v)[i]) * dx

        if geo.params["x0"]:
            Fbare0 = Fbarevol if smearMolCharge else Fbaresurf
        else: Fbare0 = None

        return Fbare0
    except:
        return None

def FbareE(geo):
    try:
        # for gradient recovery, or testing with manifactured E field
        def Fbaresurf(E, i):
            dS = geo.dS("chargeddnab")
            return dolfin.Constant(DNAqs) * (-r*E[i])('-') * lscale*dS
        def Fbarevol(E, i):
            dx = geo.dx("molecule")
            return dolfin.Constant(Moleculeqv) * (-r*E[i]) * dx

        if geo.params["x0"]:
            Fbare0 = Fbarevol if smearMolCharge else Fbaresurf
        else: Fbare0 = None

        return Fbare0
    except:
        return None


def CurrentPB(geo):
    try:
        # current approximated by linear functional of potential
        bV0 = bV if bV is not None else 0.1
        def J0(v):
            c0 = cFarad*bulkcon
            L = geo.params["l0"]/lscale
            # linear PB approximation of charge
            rho = 2*c0*v
            # assume constant external electric field
            Ez = bV0/L
            # assume current dominated by el. drift
            Jz = mu*rho*Ez # = 2*c0*bV/L * v

            return r/L*Jz*geo.dx("pore")
        return J0
    except:
        return None
