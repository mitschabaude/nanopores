import os

# directory where this file is located
INSTALLDIR = os.path.dirname(os.path.realpath(__file__))

# directory where additional data (e.g. meshes, solutions) are stored
# this expands like $HOME/.nanopores
# os.path.join is better than "/" because platform-independent
HOME = os.path.expanduser("~")
NAME = "nanopores"
DATADIR = os.path.join(HOME,".%s" %NAME)
DROPBOX = os.path.join(HOME, "Dropbox", "nanopores")
DROPBOX_FIGDIR = os.path.join(HOME, "Dropbox", "nanopores", "figures")

# this should be done by every IO function separately
#if not os.path.exists(DATADIR):
#    os.makedirs(DATADIR)

del os
