from nanopores import *

A = Box([0,0], [2,2])
B = Box([1,1], [3,3])
A_ = Box([0,1e-12], [1,1])

print A
print B

print A == A_
print A < B

print interval_union(A.intervals[0], A_.intervals[0])
print interval_union(B.intervals[0], A.intervals[0])
print interval_union(A.intervals[0], B.intervals[0])

l = BoxCollection([A_, B, A])
print l
l.sort()
print l
print
C = A + B

print C
print

D = C + A_

print D
print
print D + C
