import numpy as np

timer=np.load('timer.npy')
counter=np.load('counter.npy')
print('sims = %d'%np.sum(counter))
print('size_timer = %d'%timer.shape[0])
print(counter)
prop=float(counter[0])/float(np.sum(counter))*100
print('prop = %.2f percent' %prop)
