""" --- geometry parameters for Howorka cylindrical geometry --- """
# modified parameters to create nice pictures

nm = 1e-9

# @TODO maybe include tolc in parent file
tolc = 1e-14*nm  # tolerance for coordinate comparisons

# DNA radius
rDNA = 1.1*nm
# molecule radius
rMolecule = 0.5*nm
# effective pore radius
r0 = 2*nm #1.1*rDNA
# barrel outer radius
r1 = r0 + 2*rDNA #2.75*nm
# pore length
l0 = 15.0*nm
# membrane thickness
l1 = 4*nm
# Radius of domain
Ry = Rz = 10.0*nm
Rx = R = 10.0*nm

# total number of charged DNA base pairs
ncbp = 6.0*36  # 6 scaffold strands of DNA with 36 charged base pairs

# length scales
lc = nm
lcMolecule = lc/2
lcDNA = lc/2
lcFluid = lc
lcFluidCtr = lc/2
