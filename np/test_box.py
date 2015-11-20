from nanopores import *

A = Box([0,0,0], [2,2,2])
B = Box([1,1,0], [3,3,1])
C = Box([0.5,-0.5,0], [4,0.5,1])

ABC = Box.union(A, B, C)
ABC.plot()
