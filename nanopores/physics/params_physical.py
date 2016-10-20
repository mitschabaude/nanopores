"""
physical parameters for Simulations
"""
nm = 1e-9
cm = 1e-2
pi = 3.14159265359

# Define constants taken from
# 'Effective force applied on DNA inside a solid-state nanopore'
qq = 1.602e-19  # positive elementary charge [C]
mol = 6.022e23  # Avogadro
kB = 1.3806488e-23 # boltzmann [J/K]
T = 293 # room temperature (20 grad Celsius) [K]
eperm = 8.854e-12  # electrical permittivity [F/m]
rpermw = 80.2  # relative permittivity or dielectric constant of water
rpermSiN = 7.0  # relative permittivity of SiN (http://www.iue.tuwien.ac.at/phd/quay/node27.html)
rpermLipid = 2.0  # relative permittivity of a lipid bilayer  # @TODO
rpermDNA = 12.0  # relative permittivity of dsDNA # http://www.nature.com/nmat/journal/v11/n9/pdf/nmat3369.pdf
rpermSAM = 2.7  # according to Wei sensing
rpermAu = 6.9  # according to http://www.mit.edu/~6.777/matprops/gold.htm
D = 1.9e-9  # Diffusion [m^2/s] # lower in pore according to kurnikowa paper
eta = 1e-3  # fluid viscosity [Pa s]
rho = 1e3 # fluid density [kg/m^3]
Reynolds = rho*0.01*10*nm/eta # typical Reynolds number in nanopore systems with u = 0.01 m/s

# derived constants
kT = kB*T
UT = kB*T/qq # thermal voltage [V] = ca. 25 mV
cFarad = qq*mol  # Faraday [C/mol] cFarad=qq*mol
mu = D/UT # 73e-9  # average mobility [m^2/Vs]

bpq = -2*qq  #charge per base pair -keyser.pdf
distbp = 0.34*nm  # distance between two base pairs for dsDNA -keyser.pdf
DNAql = bpq/distbp  # line charge of dsDNA
DNAqs = -qq/nm**2 # surface charge of dsDNA (according to Howorka)

# absolute permittivities
permittivity = {
    'water':eperm*rpermw,
    'dna':eperm*rpermDNA,
    'molecule':eperm*rpermDNA,
    'sin':eperm*rpermSiN,
    'lipid':eperm*rpermLipid,
    'au':eperm*rpermAu,
    'sam':eperm*rpermSAM,
    'oxide':eperm*3.9,
    'sil':eperm*11.7,
}

