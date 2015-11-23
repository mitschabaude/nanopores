from nanopores import *

A = Box([0,0,0], [2,2,2])
B = Box([1,1,-1], [3,3,1])
C = Box([0.5,-0.5,0.5], [4,0.5,1.5])

ABC = Box.union(A, B, C)
ABC.plot()

print ABC.esets[3]

D = Box([2,1,1], [3,3,3])

ABC.addboxes(D)
ABC.plot()

print ABC.esets[3]
