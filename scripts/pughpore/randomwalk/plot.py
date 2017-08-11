import sys
if len(sys.argv)!=4:
    print '[type,traj,both], samples (int), fieldsname'
    exit()
if not (sys.argv[1]=='type' or sys.argv[1]=='traj' or sys.argv[1]=='both' or sys.argv[1]=='None'):
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
    if fields.exists(fieldsname):
        enter = raw_input("%s already exists! Continue? [y/n]}\n"%fieldsname)
        if enter != 'y'
            exit()
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
#params=dict(avgbind1=17.2e6,avgbind2=3e4,P_bind1=0.193,P_bind2=3e-1,z0=hpore/2.+0.)
params=dict(avgbind1=17.2e6,avgbind2=3e4,P_bind1=0.020,P_bind2=3e-1,z0=hpore/2.+0.)
#b1 = [[[l3/2.,-hpore/2.],[l3/2.,hpore/2.-h2],[l2/2.,hpore/2.-h2],[l2/2.,hpore/2.-h1],[l1/2.,hpore/2.-h1],[l1/2.,hpore/2.]]]
#b2 = [[[2.5, hpore/2.-h2-24.], [2.5, hpore/2.-h2-32.]]]
b2 = [[[l3/2.,-hpore/2.],[l3/2.,hpore/2.-h2],[l2/2.,hpore/2.-h2],[l2/2.,hpore/2.-h1],[l1/2.,hpore/2.-h1],[l1/2.,hpore/2.]]]
b1 = []
#b2 = []
outside=True
#outside=False

for i in range(samples):
	run(params,fieldsname,outcome,outside,b1,b2)
	print '%i out of '%i+num 
try: fields.get_fields(fieldsname,**params)["b1"]
except: fields.save_fields(fieldsname,params,b1=b1)
try: fields.get_fields(fieldsname,**params)["b2"]
except: fields.save_fields(fieldsname,params,b2=b2)

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
