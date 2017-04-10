import numpy as np
import os
import nanopores.tools.fields as f
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)

fieldsname='number_of_collisions'
params=dict(avgbind1=2e7,avgbind2=3e4,P_bind1=0.,P_bind2=0.,z0=23.)

data=f.get_fields(fieldsname,**params)
Nc=np.array(data["Nc"])
print Nc
print 'len = '+str(Nc.shape[0])
print np.mean(Nc)
