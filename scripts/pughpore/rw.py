# (c) 2017 Gregor Mitscha-Baude
"""Random walks in Pugh pore on 2D proxy geometry, to determine distribution of
attempt time needed to fit binding parameters."""
import numpy as np
import matplotlib.pyplot as plt
import dolfin
import nanopores
import nanopores.models.randomwalk as randomwalk
from nanopores.tools import fields
fields.set_dir_mega()

params = nanopores.user_params(
    # general params
    # geo
    geoname = "pughcyl",
    dim = 2,
    diamPore = 6.,
    rMolecule = 2.0779,
    R = 40.,
    Htop = 60.,
    Hbot = 35.,
    geop = dict(R=40., Htop=60., Hbot=35.),
    x0 = None,
    # physics
    Qmol = 5.,
    bulkcon = 1000.,
    dnaqsdamp = 0.7353,
    bV = -0.1,
    posDTarget = True,
    # solver
    h = 2.,
    frac = 0.5,
    Nmax = 5e4,
    imax = 30,
    tol = 1e-3,
    cheapest = False,
    stokesiter = False, #True
    hybrid = True,
    reconstruct = False,

    # random walk params
    N = 100, # number of (simultaneous) random walks
    dt = .3, # time step [ns]
    walldist = 1.2, # in multiples of radius, should be >= 1
    margtop = 35.,
    margbot = 0.,
    rstart = 5.,
    initial = "sphere",

    # receptor params
    lbind = 4., # length of binding zone for long binding [nm]
)
########### WHAT TO DO  ###########
todo = nanopores.user_params(
    test_solver = False,
    plot_dolfin = False,
    plot_streamlines = False,
    video = False,
)

########### CONSTANTS ###########
setup_rw = randomwalk.setup_default

########### TEST AND PLOT SOLVER ###########
if todo.test_solver:
    import nanopores.models.nanopore as nanopore
    setup = nanopore.Setup(**params)
    _, pnps = nanopore.solve(setup, True)    
    dolfin.interactive()

if todo.plot_dolfin:
    rw = setup_rw(params)
    dolfin.plot(rw.D[1])
    dolfin.plot(rw.F)
    dolfin.interactive()
    
if todo.plot_streamlines:
    rw = setup_rw(params)
    rw.plot_streamlines(both=True, R=20, Hbot=30, Htop=35,
                        maxvalue=1e-10, figsize=(5, 5))
    plt.figure("D")
    dolfin.plot(rw.D[1], backend="matplotlib")
    plt.show()
    
if todo.video:
    rw = setup_rw(params)
    randomwalk.run(rw)