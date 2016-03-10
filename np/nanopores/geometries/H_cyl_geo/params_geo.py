""" --- geometry parameters for Howorka cylindrical geometry --- """

nm = 1e0

# @TODO maybe include tolc in parent file
tolc = 1e-10*nm  # tolerance for coordinate comparisons

dim = 3

# DNA radius
rDNA = 1.1*nm
# molecule radius
rMolecule = 0.5*nm
# effective pore radius
r0 = 1.*nm
# barrel outer radius
r1 = 2.5*nm
# pore length
l0 = 9.0*nm
# membrane thickness
l1 = 2.2*nm
# Radius of domain
Rz = 10.0*nm
R = 10.0*nm

# total number of charged DNA base pairs
ncbp = 6.0*36  # 6 scaffold strands of DNA with 36 charged base pairs

# characteristic length / mesh size h
lc = nm
lcMolecule = lc*3e-2
lcCenter = lc*5e-2
lcOuter = lc

# provide default values for boundary layer around membrane/molecule
membraneblayer = None
moleculeblayer = None
