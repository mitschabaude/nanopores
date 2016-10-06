# (c) 2016 Gregor Mitscha-Baude
" 1D PNP, modelling reservoirs and membrane far away from pore "

import nanopores as nano
import solvers

geop = nano.Params(
    R = 35.,
    H = 70.,
)
physp = nano.Params(
    bulkcon = 1000.,
    bV = -1.,
)

geo, pnp = solvers.solve1D(geop, physp)
solvers.visualize1D(geo, pnp)
nano.showplots()

