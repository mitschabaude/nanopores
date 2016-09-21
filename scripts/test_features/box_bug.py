from nanopores import Box

solid = Box((0, 0), (3, 3)) - Box((1, 1), (2, 2))
solid.addsubdomain(solid, "solid")
#solid.addboundary(Box((0, 0), (0, 3)), "left")

#solid = (A | B) - C
print solid
geo = solid.create_geometry(lc=.2)
solid.plot()