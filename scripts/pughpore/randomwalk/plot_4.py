import sys
if len(sys.argv)!=4:
    print '[type,traj,both], samples (int), fieldsname'
    exit()
if not (sys.argv[1]=='type' or sys.argv[1]=='traj' or sys.argv[1]=='both'):
    print '[type,traj,both], samples (int), fieldsname'
    exit()
else:
    outcome=sys.argv[1]
try:
    samples=int(sys.argv[2])
    num = sys.argv[2]
except:
    print '[type,traj,both], samples (int), fieldsname'
    exit()
fieldsname = sys.argv[3]
if samples!=0:
	from run import *
else:
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
#params=dict(avgbind1=7e5,avgbind2=1e3,P_bind1=1.e-2,P_bind2=1e-1,z0=hpore/2.+5.)
#params=dict(avgbind1=1e7,avgbind2=3e4,P_bind1=5.e-3,P_bind2=8e-2,z0=hpore/2.+0.)
#params=dict(avgbind1=1e7,avgbind2=1e5,P_bind1=2.e-2,P_bind2=8e-2,z0=hpore/2.+0.) # good if only real translocations count - no type 0 "ood"
params=dict(avgbind1=2e7,avgbind2=3e4,P_bind1=8.e-2,P_bind2=0*3e-1,z0=hpore/2.+0.)
b1 = []
b2 = [[[l3/2.,-23.],[l3/2.,-9.]]]
outside=True

for i in range(samples):
	run(params,fieldsname,outcome,outside,b1,b2)
	print '%i out of '%i+num 

print 'field updates'
if outcome=='type' or outcome=='both':
    from create_plot_type import *
    f.update()
    save_fig_type(params,fieldsname)
    plt.close("all")
    from create_plot_traj import *
    save_fig_traj(params,fieldsname,0,False)

if outcome=='traj' or outcome=='both':
    if outcome!='both':
        from create_plot_traj import *
        f.update()
    for i in range(len(fields.get_fields(fieldsname,**params)["X"])):
        save_fig_traj(params,fieldsname,i,True)
