import numpy as np
file=open('taus','r')
for line in file:
    tau=line.split('\n')[0]
    data = np.load('data'+str(tau)+'.npy')
    print tau
    print np.mean(data)
