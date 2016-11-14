# (c) 2016 Gregor Mitscha-Baude
import os

HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")

DATADIR = os.path.join(HOME, "Dropbox", "Paper Howorka", "data", "fields")

import nanopores.tools.fields as fields
fields.set_dir(DATADIR)
