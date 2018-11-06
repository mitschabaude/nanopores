# (c) 2016 Gregor Mitscha-Baude
import os

HOME = os.path.expanduser("~")
#PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
#FIGDIR = os.path.join(PAPERDIR, "figures", "")

FIGDIR = os.path.join(HOME, "Dropbox", "nanopores", "figures")
FIGDIR_HOWORKA = os.path.join(HOME, "Dropbox", "Paper Howorka", "figures")
#DATADIR = os.path.join(HOME, "Dropbox", "nanopores", "fields")

import nanopores.tools.fields as fields
fields.set_dir_mega()
#fields.set_dir(DATADIR)
