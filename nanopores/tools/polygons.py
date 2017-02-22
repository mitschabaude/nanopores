# (c) 2017 Gregor Mitscha-Baude
"""this was intended to perform boolean operations on polygons,
but this seems to be too hard without relying on external software"""

def clip(subjectPolygon, clipPolygon):
   def inside(p):
      return(cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0])

   def computeIntersection():
      dc = [ cp1[0] - cp2[0], cp1[1] - cp2[1] ]
      dp = [ s[0] - e[0], s[1] - e[1] ]
      n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
      n2 = s[0] * e[1] - s[1] * e[0]
      n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0])
      return [(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3]

   outputList = subjectPolygon
   cp1 = clipPolygon[-1]

   for clipVertex in clipPolygon:
      cp2 = clipVertex
      inputList = outputList
      outputList = []
      s = inputList[-1]

      for subjectVertex in inputList:
         e = subjectVertex
         if inside(e):
            if not inside(s):
               outputList.append(computeIntersection())
            outputList.append(e)
         elif inside(s):
            outputList.append(computeIntersection())
         s = e
      cp1 = cp2
   return(outputList)

def horizontal_intersect(poly, a, b):
    poly2 = [[-100, a], [-100, b], [100, b], [100, a]]
    return clip(poly2, poly)

if __name__ == "__main__":
    poly1 = [[0,0], [2,0], [2,2], [0,2]]
    poly2 = [[1,1], [3,1], [3,3], [1,3]]
    poly3 = [[0.5,-1], [0.5,3], [1.5,3], [1.5,-1]]
    u = [[0,-2], [0,2], [0.5,2], [0.5,-1], [1.5,-1], [1.5,2], [2,2], [2,-2]]

    #print "p1", poly1
    #print "p2", poly2
    #print "intersection", clip(poly1, poly2)
    #
    #print clip(u, poly1)

    print horizontal_intersect(poly1, .5, 1.5)