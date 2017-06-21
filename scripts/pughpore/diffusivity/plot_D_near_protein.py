# (c) 2017 Gregor Mitscha-Baude
import dolfin
import os, sys
import numpy as np
from nanopores.models.pughpore import *
from nanopores import user_params

setup = Setup(**user_params(dim=3, h=2., Nmax=5e5, x0=[0.,0.,0.],
               diffusivity="Dpugh2", cheapest=True, diamPore=6.))
setup.prerefine()
set_D(setup)
dim = setup.phys.dim
setup.plot(setup.phys.Dp[dim-1, dim-1], interactive=True)
