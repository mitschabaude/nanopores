""" simple base class for pore geometries compatible with PNPS """

import dolfin
from nanopores.physics.electrolyte import *

bV = None # voltage bias across pore [V] (None for no enforced bias)

Membraneqs = -0.0 # membrane surface charge density [C/m**2]
DNAqsPure = -qq/nm**2 # = -0.16 # DNA surface charge density [C/m**2]
dnaqsdamp = 1. # DNA charge damping
SiNqs = -0.022
SAMqs = -0.078
ahemqs = 0.

rpermPore = rpermw
rpermProtein = 2. # TODO ?????
rDPore = 0.5

DNAqs = lambda: DNAqsPure*dnaqsdamp
permPore = lambda: eperm*rpermPore
permProtein = lambda: eperm*rpermProtein
DPore = lambda: D*rDPore

# piece-wise boundary conditions
v0 = dict(
    upperb = 0.,
    lowerb = "bV"
)
c0 = dict(
    upperb = "bulkcon",
    lowerb = "bulkcon"
)
cp0 = cm0 = c0

permittivity.update(
    #default = eperm*rpermw,
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
    #default = "D",
    bulkfluid = "D",
    pore = "DPore",
    solid = 0.,
)
Dm = Dp

# TODO: this seems a little much; clean it up
synonymes = {
    "pore": {"poretop", "porecenter", "porebottom"},
    "bulkfluid": {"fluid_bulk_top", "fluid_bulk_bottom"},
    "fluid": {"pore", "bulkfluid"},
    "solid": {"membrane", "channel", "molecule"},
    "protein": {"ahem"},
    "channel": {"ahem"},
    "proteinb": {"ahemb"},
    "noslip": {"proteinb", "membraneb"},
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

