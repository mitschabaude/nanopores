import numpy as np
import os
import nanopores.tools.fields as fields
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
fields.set_dir(DATADIR)
old=['events3_onlyone_%i' %i for i in [1,2,3,4,5]]
new=['events3_onlyone_%i_new' %i for i in [1,2,3,4,5]]

hpore=46.
params=dict(avgbind1=2e7,avgbind2=3e4,P_bind1=8.e-2,P_bind2=0*3e-1,z0=hpore/2.+0.)
for i in range(5):
    print(i)
    data=fields.get_fields(old[i],**params)
    for k in range(len(data["X"])):
        if k%10==0:
            print(k)
        saveX=[np.array(data["X"][k])]
        saveY=[np.array(data["Y"][k])]
        saveZ=[np.array(data["Z"][k])]
        fields.save_fields(new[i],params,X=saveX,Y=saveY,Z=saveZ)
    fields.save_fields(new[i],params,a=data["a"],ood=data["ood"],J=data["J"],t=data["t"],b1=data["b1"],Fzavg=data["Fzavg"],b2=data["b2"],Dzavg=data["Dzavg"],T=data["T"])
    fields.update()
