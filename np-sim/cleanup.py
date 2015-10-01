#!/usr/bin/env python

''' clean up data file in case of erroneous calculation '''

from nanopores import DATADIR
from protocol import Data

# get git repo
with open("%s/sim/gitrepo"%(DATADIR,), "r") as f:
    gitrepo = f.next().strip()
    
# get data filename
with open("%s/sim/datafile"%(DATADIR,), "r") as f:
    datafile = f.next().strip()
    
folder = "np-sim/simulations"
    
# access data file
FILENAME = "/".join([gitrepo, folder, datafile])
data = Data(FILENAME)

# CLEAN UP: set free all vertices that are currently in calculation
data.set_free()




