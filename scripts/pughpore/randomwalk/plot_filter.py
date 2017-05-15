import os
import sys
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)
start = int(sys.argv[2])
fieldsname = sys.argv[1]

params=dict(avgbind1=23e6,avgbind2=3e4,P_bind1=0.035,P_bind2=3e-1,z0=hpore/2.+0.)
data=f.get_fields(fieldsname,**params)
samples=len(data["a"])
from create_plot_filter import save_fig_filter
for i in range(start,samples):
    print '%i out of %i' %(i,samples)
    save_fig_filter(params,fieldsname,i)

f.update()
