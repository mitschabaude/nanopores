from nanopores import *

A = Box([0,0,0], [2,2,2])
B = Box([1,1,0], [3,3,1])
C = Box([0.5,-0.5,0], [4,0.5,1])
A_ = Box([0,1e-12], [2,2])

#print A == A_ # True
#print A < B # True

#AB = Box.union(A, B)
#print AB
#AC = Box.union(A, C)
#print AC

ABC = Box.union(A, B, C)

#print ABC
ABC.plot()
