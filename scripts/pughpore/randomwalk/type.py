import sys
num = sys.argv[1]
print 'run RW'
if int(num)!=0:
	from new_run import *
else:
	import nanopores as nano
	import nanopores.geometries.pughpore as pughpore
	geop = nano.Params(pughpore.params)
	hpore = geop.hpore
#params=dict(avgbind1=7e5,avgbind2=1e3,P_bind1=1.e-2,P_bind2=1e-1,z0=hpore/2.+5.)
params=dict(avgbind1=1e7,avgbind2=3e4,P_bind1=5.e-3,P_bind2=8e-2,z0=hpore/2.+5.)

for i in range(int(num)):
	run(params)
	print '%i out of %i'%(i,int(num))

print 'field updates'
from plot3 import *
f.update()
save_fig(params)
