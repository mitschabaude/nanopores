# -*- coding: utf-8 -*-
import nanopores as nano
import numpy as np
import os
import nanopores.geometries.pughpore as pughpore
import nanopores.tools.fields as f
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)

number=False

geop = nano.Params(pughpore.params)
hpore=geop.hpore
fieldsname='events10_onlyone'
params=dict(avgbind1=2e7,avgbind2=3e4,P_bind1=1.074/1.74,P_bind2=0*3e-1,z0=hpore/2.+0.)
data=f.get_fields(fieldsname,**params)
t = np.array(data["t"])
a= np.where(t>0.2)[0].shape[0]
b= t.shape[0]
print a
print b
print float(a)/float(b)
