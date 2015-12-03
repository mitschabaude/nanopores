import numpy as np

counter_1 = np.load('counter_1.npy')
counter_2 = np.load('counter_2.npy')
counter_3 = np.load('counter_3.npy')
counter_4 = np.load('counter_4.npy')
counter_5 = np.load('counter_5.npy')
counter_6 = np.load('counter_6.npy')

counter=np.zeros_like(counter_1)

for index in range(counter.shape[0]):
	counter[index] = np.sum(np.array([counter_1[index],counter_2[index],counter_3[index],counter_4[index],counter_5[index],counter_6[index]]))



exit_x = np.concatenate((np.load('exit_x_1.npy'),np.load('exit_x_2.npy'),np.load('exit_x_3.npy'),np.load('exit_x_4.npy'),np.load('exit_x_5.npy'),np.load('exit_x_6.npy')))
exit_y = np.concatenate((np.load('exit_y_1.npy'),np.load('exit_y_2.npy'),np.load('exit_y_3.npy'),np.load('exit_y_4.npy'),np.load('exit_y_5.npy'),np.load('exit_y_6.npy')))
exit_z = np.concatenate((np.load('exit_z_1.npy'),np.load('exit_z_2.npy'),np.load('exit_z_3.npy'),np.load('exit_z_4.npy'),np.load('exit_z_5.npy'),np.load('exit_z_6.npy')))

np.save('exit_x',exit_x)
np.save('exit_y',exit_y)
np.save('exit_z',exit_z)
np.save('counter',counter)

