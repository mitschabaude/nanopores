# (c) 2016 Gregor Mitscha-Baude
import numpy as np
from nanopores.geometries.pughpore import params as pugh_params
import nanopores.models.pughpore as pugh
from nanopores import Params


setup = pugh.SetupNoGeo()
params = Params(pugh_params) | setup.geop
#print params

# ---- create z part of tensor grid -----
r = params.rMolecule
eps = 1e-2
r = r + eps
buf = 5.

ztop = params.hpore/2.
zbot = -ztop

# 6 nodes => 5 sections
z = [zbot - buf,
     zbot - r,
     ztop - params.h2 + r,
     ztop - params.h1 + r,
     ztop + r,
     ztop + buf]

# lengths of sections, relative meshwidth, radii
lengths = np.diff(np.array(z))
hz = np.array([1.5, 1, 1, 1, 1.5])
rpore = [params.l0/2.,
         params.l3/2. - r,
         params.l2/2. - r,
         params.l1/2. - r,
         params.l0/2.]


print z
print lengths



#........................R.............................
#                                                     .
#                                                     .
#              .........l0..........                  .
#              .                   .                  .
#              ._ _______________ _...............    .
#              |D|               |D|     .   .   .    .
#              |D|......l1.......|D|    h1   .   .    .
#              |D|_ ____l2_____ _|D|......   h2  .    .
#              |DDD|_ _______ _|DDD|..........   .    .
#              |DDDDD|       |DDDDD|             .    .
#              |DDDDD|       |DDDDD|             .    .
#       DNA--->|DDDDD|       |DDDDD|           hpore  .
#              |DDDDD|       |DDDDD|             .    .
#              |DDDDD|..l3...|DDDDD|             .    .
#   MEMBRANE   |DDDDD|       |DDDDD|             .    H
#      |       |DDDDD|       |DDDDD|             .    .
#      |       |DDDDD|       |DDDDD|....h4       .    .
#______V_________|DDD|       |DDD|_____.________ .___ .......
#MMMMMMMMMMMMMMMM|DDD|       |DDD|MMMMM.MMMMMMMMM.MMMM.    hmem
#MMMMMMMMMMMMMMMM|DDD|_______|DDD|MMMMM.MMMMMMMMM.MMMM.......
#                .               .                    .
#                .......l4........                    .
#                                                     .
#                                                     .
#                                                     .
#......................................................