from nanopores import Box

#zero = [0, 0, 0]
#
#R = 30
#H = 40
#hmem = 2
#rdna = 20
#hdna = 20
#cdna = [0, 0, 0.5*(hdna-hmem)]
#rpore = 10
#
#A = Box(center=zero, l=R, w=R, h=hmem)
#B = Box(center=cdna, l=rdna, w=rdna, h=hdna)
#C = Box(center=cdna, l=rpore, w=rpore, h=hdna)

solid = Box((0, 0), (3, 3)) - Box((1, 0), (2, 2))
#solid.addsubdomain(solid, "solid")
#solid.addboundary(Box((0, 0), (0, 3)), "left")

#solid = (A | B) - C
print solid
geo = solid.create_geometry(lc=.2)
solid.plot()