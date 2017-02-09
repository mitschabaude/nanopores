# (c) 2016 Gregor Mitscha-Baude
import os

HOME = os.path.expanduser("~")
#PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
#FIGDIR = os.path.join(PAPERDIR, "figures", "")

FIGDIR = os.path.join(HOME, "Dropbox", "nanopores", "figures")
DATADIR = os.path.join(HOME, "Dropbox", "nanopores", "fields")

import nanopores.tools.fields as fields
fields.set_dir(DATADIR)
