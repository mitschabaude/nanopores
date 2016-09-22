from nanopores.tools.box import Box

square = Box((0, 0), (4, 4))
lshape = Box((1, 1), (3, 3)) - Box((2, 2), (3, 3))

square.addsubdomain(lshape, "lshape")
square.addboundaries(   
    allb = (square - lshape).boundary(),
)
geo = square.create_geometry(lc=.3)

square2 = square.copy()
geo = square2.create_geometry(lc=.3, merge=False)

square.plot()
square2.plot()