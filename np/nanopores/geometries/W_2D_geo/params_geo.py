""" --- geometry parameters for Wei et al stochastic sensing --- """

nm = 1e-9

# @TODO maybe include tolc in parent file
tolc = 0.01*nm  # tolerance for coordinate comparisons

dim = 2

# molecule radius
rMolecule = 5*nm
# effective pore radius ~ d^SAM_p/2
r0 = 13*nm

# aperture angle in degrees
angle = 40

# SiN membrane thickness (in vertical direction)
lsin = 50*nm
# Au membrane thickness (in vertical direction)
lau = 40*nm
# Au thickness in radial direction
rlau = 10*nm
# SAM layer thickness (in vertical direction)
lsam = 3*nm

# Radius of domain
Rz = 150.0*nm
R = 150.0*nm

# fraction of R which is considered as outer membrane boundary
outerfrac = 0.2

# mesh generation parameters
# length scales
lc = 10*nm
lcMolecule = lc*1e-1
lcOuter = lc*2
lcCenter = lc/5

boxfields = None
# provide default values for boundary layer around membrane/molecule
moleculeblayer = None
# if True a boundary layer field is applied, elsewise a box field
membraneblayer = None
lowerblayer = None
