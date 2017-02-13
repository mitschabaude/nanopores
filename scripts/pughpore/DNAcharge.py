# (c) 2017 Gregor Mitscha-Baude
from math import pi, sqrt

charge = -2 # charge per base pair
distbp = 0.34 # distance between two base pairs for dsDNA
diam = 2.5 # DNA diameter

Abox = diam * 4 * distbp # surface area of DNA box
Acyl = diam * pi * distbp # surface area of DNA cylinder
Acylinner = 3 * diam*2./sqrt(pi) * pi * distbp

DNAqsbox = charge / Abox
DNAqscyl = charge / Acyl
DNAqsinner = 3 * charge / Acylinner

print "DNA charge (box): %.4f [q/nm^2]" % DNAqsbox
print "DNA charge (cyl): %.4f [q/nm^2]" % DNAqscyl
print "DNA charge (inner cyl): %.4f [q/nm^2]" % DNAqsinner
