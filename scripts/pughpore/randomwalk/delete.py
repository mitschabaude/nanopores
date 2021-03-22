# (c) 2016 Gregor Mitscha-Baude
import os
import sys

HOME = os.path.expanduser("~")
#PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
#FIGDIR = os.path.join(PAPERDIR, "figures", "")

FIGDIR = os.path.join(HOME, "Dropbox", "nanopores", "figures")
DATADIR = os.path.join(HOME, "Dropbox", "nanopores", "fields")

import nanopores.tools.fields as fields
fields.set_dir(DATADIR)
fieldsname=sys.argv[1]
if fields.exists(fieldsname):
	a=input('delete %s'%fieldsname+' ?')
	if a=='y' or a=='Y':
		fields.remove(fieldsname)
