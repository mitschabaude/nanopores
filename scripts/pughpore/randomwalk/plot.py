import sys
if len(sys.argv)!=4:
    print '[type,traj], samples (int), fieldsname'
    exit()
if not (sys.argv[1]=='type' or sys.argv[1]=='traj'):
    print '[type,traj], samples (int), fieldsname'
    exit()
else:
    outcome=sys.argv[1]
try:
    samples=int(sys.argv[2])
    num = sys.argv[2]
except:
    print '[type,traj], samples (int), fieldsname'
    exit()
#try:
#    fieldsname = sys.argv[3]
#except:
#    fieldsname = 'rw_type_test'
if samples!=0:
	from run_type import *
else:
	import nanopores as nano
	import nanopores.geometries.pughpore as pughpore
	geop = nano.Params(pughpore.params)
	hpore = geop.hpore
#params=dict(avgbind1=7e5,avgbind2=1e3,P_bind1=1.e-2,P_bind2=1e-1,z0=hpore/2.+5.)
params=dict(avgbind1=1e7,avgbind2=3e4,P_bind1=5.e-3,P_bind2=8e-2,z0=hpore/2.+0.)

for i in range(samples):
	run(params,fieldsname,outcome)
	print '%i out of '%i+num 

print 'field updates'
if outcome=='type':
    from create_plot_type import *
    f.update()
    save_fig(params,fieldsname)
elif outcome=='traj':
    from create_plot_traj import *
    f.update()
    save_fig(params,fieldsname,0)
