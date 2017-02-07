from run import *
params=dict(avgbind=1e7,P_bind=3.e-4,z0=hpore/2.+5.)
num = sys.argv[1]

for i in range(int(num)):
	run(params)
