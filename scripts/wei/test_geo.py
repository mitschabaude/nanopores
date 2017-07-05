# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from nanopores.tools.polygons import Ball, Polygon, MultiPolygon, MultiPolygonPore
from nanopores.geometries.cylpore import MultiPore
from nanopores import user_params
params = user_params(
    R = 100.,
    R0 = 100.,
    H = 150.,
    cs = [-15.5, 15.5],
    x0 = [0, 0, 10],
    rMolecule = 2.1,
    dim = 3,
    no_membrane = True,
    r0 = 13, # pore radius
    angle = 40, # aperture angle in degrees
    lcCenter = 0.3,
    lcMolecule = 0.1,

)
# SiN membrane thickness (in vertical direction)
lsin = 50.
# Au membrane thickness (in vertical direction)
lau = 40.
# Au thickness in radial direction
rlau = 10.
# SAM layer thickness (in vertical direction)
lsam = 3

l0 = lau + lsin + lsam
angle2 = params.angle/2. * np.pi/180.
tan = np.tan(angle2)
sin = np.sin(angle2)
cos = np.cos(angle2)

l = l0/2.
r0 = params.r0
r1 = r0 + l0*tan
rsam = r0 + lsam/cos
rsin = r0 + lsam/cos + rlau
R = params.R0

sam = [[r0, -l], [r1, l], [R, l], [R, l - lsam],
       [rsam - tan*(lsam - l0), l - lsam], [rsam, -l]]
au = [sam[5], sam[4], sam[3], [R, -l + lsin], [rsin + tan*lsin, -l + lsin],
      [rsin, -l]]
sin = [au[5], au[4], au[3], [R, -l]]

p = MultiPore(**params)
p.add_polygons(sam=sam, au=au, sin=sin)

receptor = Ball([0.,0.,25.], 4., 2, lc=0.25)
p.add_balls(receptor=receptor)
geo = p.build_solid(h=10.)
P = p.protein
P.plot(".k")
from matplotlib import pyplot as plt
plt.xlim(0, R + 5)

print geo

geo.plot_subdomains()
geo.plot_boundaries(interactive=True)