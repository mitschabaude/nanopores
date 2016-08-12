""" base class for all electrolyte related physics
i.e. the most simple specifications that pnps should work with """

import dolfin
from nanopores.physics.default import *

T = 293 # temperature [K]
bulkcon = 300. # bulk concentration of ions [mol/m**3]
D = 1.9e-9  # diffusivity [m^2/s]

kT = lambda T: kB*T
UT = lambda kT: kT/qq
mu = lambda D, kT: D*qq/kT # mobility [m^2/Vs]
cFarad = qq*mol  # Faraday constant [C/mol]
debye = lambda bulkcon, kT: dolfin.sqrt(rpermw*eperm*kT/qq**2/2/mol/bulkcon)  # debye length [m]
bulkconduct = lambda bulkcon, mu: 2.*bulkcon*cFarad*mu # electrolyte bulk conductivity [S/m]


# piece-wise boundary conditions
v0 = dict()
c0 = dict()
cp0 = cm0 = c0

# rhs data
surfcharge = dict() # surface charge densities for Neumann RHS
volcharge = dict(default = 0.) # volume charges for RHS
cpflux = dict()
cmflux = dict()

Dp = Dm = dict(default = "D", solid = 0.)

# no-slip velocity
U0 = lambda dim: dolfin.Constant(tuple(0. for i in range(dim)))
noslip = dict(noslip = "U0")
pressure = dict(nopressure = 0.)

synonymes = dict(
    nocbc = set(),
)
