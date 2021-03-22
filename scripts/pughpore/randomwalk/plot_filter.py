import os
import sys
start = int(sys.argv[2])
fieldsname = sys.argv[1]
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

params=dict(avgbind1=23e6,avgbind2=3e4,P_bind1=0.035,P_bind2=3e-1,z0=hpore/2.+0.)
data=fields.get_fields(fieldsname,**params)
samples=len(data["a"])
from create_plot_filter import save_fig_filter
#for i in [start]:
for i in range(start,samples):
    print('%i out of %i' %(i,samples))
    save_fig_filter(params,fieldsname,i)

fields.update()
