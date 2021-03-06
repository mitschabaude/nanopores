""" --- geometry parameters for Howorka --- """

nm = 1e-9

# @TODO maybe include tolc in parent file
tolc = 1e-14*nm  # tolerance for coordinate comparisons

dim = 2

# DNA radius
rDNA = 1.1*nm
# molecule radius
rMolecule = 0.5*nm
# effective pore radius
r0 = 1*nm
# barrel outer radius
r1 = 2.5*nm
# pore length
l0 = 15.0*nm
# membrane thickness
l1 = 2.2*nm

# Radius of domain
Rz = 15.0*nm
R = 15.0*nm
Rx = R
Ry = Rz

# total number of charged DNA base pairs
ncbp = 6.0*36  # 6 scaffold strands of DNA with 36 charged base pairs

# length scales
lc = nm
lcMolecule = lc*0.1
lcCenter = lc/2
lcOuter = lc

boxfields = True
# provide default values for boundary layer around membrane/molecule
membraneblayer = None
moleculeblayer = None
