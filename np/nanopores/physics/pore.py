""" simple base class for pore geometries compatible with PNPS """

import dolfin
from nanopores.physics.default import *

bV = None # voltage bias across pore [V] (None for no enforced bias)

Membraneqs = -0.0 # membrane surface charge density [C/m**2]
DNAqsPure = -qq/nm**2 # = -0.16 # DNA surface charge density [C/m**2]
dnaqsdamp = 1. # DNA charge damping
SiNqs = -0.022
SAMqs = -0.078
ahemqs = 0.

T = 293 # temperature [K]
bulkcon = 300. # bulk concentration of ions [mol/m**3]
D = 1.9e-9  # diffusivity [m^2/s]

rpermPore = rpermw
rpermProtein = 2. # TODO ?????
rDPore = 0.5

DNAqs = lambda: DNAqsPure*dnaqsdamp
permPore = lambda: eperm*rpermPore
permProtein = lambda: eperm*rpermProtein
kT = lambda: kB*T
UT = lambda: kB*T/qq
mu = lambda: D*qq/(kB*T) # mobility [m^2/Vs]
cFarad = lambda: qq*mol  # Faraday constant [C/mol]
debye = lambda: dolfin.sqrt(rpermw*eperm*kB*T/qq**2/2/mol/bulkcon)  # debye length [m]
bulkconduct = lambda: 2.*bulkcon*qq*mol*D*qq/(kB*T)  # 2*c0*cFarad*mu # electrolyte bulk conductivity [S/m]

# piece-wise boundary conditions
v0 = dict(
    upperb = 0.,
    lowerb = "bV"
)
c0 = dict(
    upperb = "bulkcon",
    lowerb = "bulkcon"
)

permittivity.update(
    default = eperm*rpermw,
    bulkfluid = eperm*rpermw,
    pore = "permPore",
    protein = "permProtein",
)

surfcharge = dict( # surface charge densities for Neumann RHS
    chargedmembraneb = "Membraneqs",
    chargeddnab = "DNAqs",
    chargedsinb = "SiNqs",
    chargedsamb = "SAMqs",
    ahemb = "ahemqs",
)

volcharge = dict( # volume charges for RHS
    default = 0.,
)

Dp = dict(
    default = "D",
    bulkfluid = "D",
    pore = rDPore*D,
    solid = 0.,
    lipid = 0.
)

Dm = Dp

synonymes = {
    "pore": {"poretop", "porecenter", "porebottom"},
    "bulkfluid": {"fluid_bulk_top", "fluid_bulk_bottom"},
    "fluid": {"pore", "bulkfluid"},
    "solid": {"membrane", "channel", "molecule"},
    "protein": {"ahem"},
    "channel": {"ahem"},
    "proteinb": {"ahemb"},
    "noslip": {"ahemb", "membraneb", "moleculeb"},
    "bulk": {"upperb", "lowerb"},
    "nopressure": {"bulk"},
    "ground": {"upperb"},
    "bV": {"lowerb"},
    "ions": {"fluid"},
    "lipid": {"membrane"},
    "exittime": {"fluid"},
    "exit": {"poreexit"},
    "sideb": {"uppersideb", "lowersideb"},
    "upperbulkb": {"upperb", "uppersideb"},
    "lowerbulkb": {"lowerb", "lowersideb"},
}

