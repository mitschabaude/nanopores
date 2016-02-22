""" base class for all electrolyte related physics
i.e. the most simple specifications that pnps should work with """

import dolfin
from nanopores.physics.default import *

T = 293 # temperature [K]
bulkcon = 300. # bulk concentration of ions [mol/m**3]
D = 1.9e-9  # diffusivity [m^2/s]

kT = lambda: kB*T
UT = lambda: kB*T/qq
mu = lambda: D*qq/(kB*T) # mobility [m^2/Vs]
cFarad = lambda: qq*mol  # Faraday constant [C/mol]
debye = lambda: dolfin.sqrt(rpermw*eperm*kB*T/qq**2/2/mol/bulkcon)  # debye length [m]
bulkconduct = lambda: 2.*bulkcon*qq*mol*D*qq/(kB*T)  # 2*c0*cFarad*mu # electrolyte bulk conductivity [S/m]

# piece-wise boundary conditions
v0 = dict()
c0 = dict()
cp0 = cm0 = c0

surfcharge = dict() # surface charge densities for Neumann RHS
volcharge = dict( # volume charges for RHS
    default = 0.,
)

Dp = dict(
    default = "D",
    solid = 0.,
)

Dm = Dp

