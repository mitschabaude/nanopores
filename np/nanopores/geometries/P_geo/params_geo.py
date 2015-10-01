"""
parameters for a geometry that is inspired by the PRE article:
'Effective force applied on DNA inside a solid-state nanopore'
"""

nm = 1e-9

# tolerance for coordinate comparisons
tolc = 0.01*nm

dim = 2

# pore radius min @center, max
r0=5*nm
r1=20*nm

# pore length @center, membrane thickness
l0=20*nm
l1=50*nm

# Radius of Omega
Ry = 100*nm
Rz = Ry  # bad style, here only 'hard' constants should be allowed
Rx = 25*nm
R = Rx # bad style, here only 'hard' constants should be allowed

#DNA size
DNAlength=100*nm  # one Kuhn-length
DNAradius=1.1*nm
