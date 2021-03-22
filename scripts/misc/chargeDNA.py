# (c) 2016 Gregor Mitscha-Baude
"surface charge of DNA"

nm = 1.
q = 1. # positive elementary charge [q] = 1.602e-19 [C]
pi = 3.141592

bpq = -2.*q  # DNA charge per base pair
distbp = 0.34*nm  # distance between two base pairs for dsDNA
rDNA = 1.25*nm # dsDNA radius

ql = bpq/distbp # line charge

# surface charge [q/nm**2]
qsCircle = ql/(2.*pi*rDNA) # circular model
qsSquare = ql/(8.*rDNA) # square model

print(qsCircle, qsSquare)
