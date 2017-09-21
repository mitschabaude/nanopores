import nanopores as nano
import nanopores.geometries.pughpore as pughpore
geop = nano.Params(pughpore.params)
hpore = geop.hpore
import os
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")

DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")

import nanopores.tools.fields as fields
fields.set_dir(DATADIR)
l0 =        geop.l0
l1 =        geop.l1
l2 =        geop.l2
l3 =        geop.l3
l4 =        geop.l4
hpore =     geop.hpore
hmem =      geop.hmem
h2 =        geop.h2
h1 =        geop.h1
h4 =        geop.h4
rMolecule = geop.rMolecule
params=dict(avgbind1=17.2e6,avgbind2=3e4,P_bind1=0.193,P_bind2=3e-1,z0=hpore/2.+0.)
fieldsname='eventspara_onlyone_2'

data=fields.get_fields(fieldsname,**params)
#t = data["t"]
#a = data["a"]
