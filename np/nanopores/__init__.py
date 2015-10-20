from physics import *
from tools import *
from geo2xml import generate_mesh

import os

# directory where this file is located
INSTALLDIR = os.path.dirname(os.path.realpath(__file__))

# directory where additional data (e.g. meshes, solutions) are stored
# this expands like $HOME/.nanopores
# os.path.join is better than "/" because platform-independent

HOME = os.path.expanduser("~")
NAME = "nanopores"
DATADIR = os.path.join(HOME,".%s" %NAME)

# not necessary:
#if not os.path.exists(DATADIR):
#    os.makedirs(DATADIR)

from scripts import simulation2D
