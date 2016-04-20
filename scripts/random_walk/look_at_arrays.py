import numpy as np

exit_x=np.load('exit_x.npy')
exit_y=np.load('exit_y.npy')
exit_z=np.load('exit_z.npy')
timer=np.load('timer.npy')
counter=np.load('counter.npy')
print 'size_x = %d'%exit_x.shape[0]
print 'size_y = %d'%exit_y.shape[0]
print 'size_z = %d'%exit_z.shape[0]
print 'size_timer = %d'%timer.shape[0]
print counter
