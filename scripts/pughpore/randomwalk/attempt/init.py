import numpy as np
file=open('taus','r')
for line in file:
    tau=line.split('\n')[0]
    np.save('data2'+tau,np.array([]))

