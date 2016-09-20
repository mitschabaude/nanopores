"2D/3D boxes plus balls"

"""Interface:
A = Box([0, 0], [1, 1])
m = [0, 0]; r = 0.25
B = Ball(m, r)
C = A | B
C.create_geometry(lc=0.1)
"""

import box

class Ball(BoxCollection):
    pass
    
if __name__ == "__main__":
    # unit square
    A = Box([0, 0], [1, 1])
    # circle with r=0.5 centered at origin
    B = Ball([0, 0], 0.5)
    # union
    C = A | B
    
    C.create_geometry(lc=0.1)
    C.plot()