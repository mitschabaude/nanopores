import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

z=np.load('x.npy')
y=np.zeros(z.shape[0])-5.41
zeros=np.zeros(z.shape[0])
x=np.arange(z.shape[0])
plt.plot(x,z)
plt.plot(x,y,color='black')
plt.scatter(x,z)
plt.show()
