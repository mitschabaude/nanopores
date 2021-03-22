# (c) 2018 Gregor Mitscha-Baude
"""back-of-the-envelope calculation of continuous vs discrete ion concnetration
at given voltage, when ion radius is chosen consistently with Stokes' law

TODO: eigentlich müsste man das in einer realistischen geometrie ausrechnen,
wo nur die surface charge vorgegeben werden muss, und dann die durchschn.
netto-ladungsdichte (c+ - c-) 1, 3 und 5 ionenradii entfernt berechnen
"""

from math import pi, exp

# physical constants, in SI
k = 1.38064e-23
T = 273 + 20
kT = k*T
mol = 6.022e+23
eta = 1e-3
q = 1.602e-19

# we use length units of nm, time units of ns => diffusivity in nm**2/ns
nm = 1e-9
ns = 1e-9

# quantities derived from stokes' law
def Dstokes(r):
    r *= nm
    return kT/(6*pi*eta*r) / (nm**2/ns)
def rstokes(D):
    D *= nm**2/ns
    return kT/(6*pi*eta*D) / nm

# experimental data
# name : measured radius [nm], diffusivity [nm**2/ns], charge valency
data = {
    "K+" : dict(r=0.152, D=1.96, z=1),
    "Na+" : dict(r=0.116, D=1.33, z=1),
    "Cl-" : dict(r=0.167, D=2.03, z=-1)
}

# volume in nm**3
def volume(r):
    return 4.*pi/3 * r**3

# volume of ions in nm**3, consistent with stokes law
vol = {ion : volume(rstokes(data[ion]["D"])) for ion in data}

# millimole and mole, in numbers per nm**3
mM = mol / (1/nm**3)
M = 1e3 * mM

# positive ion concentration (in 1/nm**3), according to boltzmann distribution,
# as a function of ion type, potential (in V), bulk concentration (in M)
def c(ion="K+", phi=-0.1, c0=0.15):
    Q = q * data[ion]["z"]
    c0 *= M    
    # boltzmann
    return c0 * exp(-Q*phi/kT)
    
# numbers of ions in an ion
def number(ion="K+", phi=-0.1, c0=0.15):
    return vol[ion] * c(ion, phi, c0)

# net charge (in q) for symmetric KCl in ion
def net_charge(phi=-0.1, c0=0.15):
    return number("K+", phi, c0) - number("Cl-", phi, c0)

print("charge [c] in one ion in 1M KCl at -50mV:", net_charge(-0.05, 1))
