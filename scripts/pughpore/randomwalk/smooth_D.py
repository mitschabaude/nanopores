from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
D=np.load('D.npy')
Z=np.load('D_Z.npy')

a, b = -7.3, 7.3
D_N = np.zeros(Z.shape[0])
for i in range(Z.shape[0]):
    if Z[i] >a and Z[i]<b:
        D_N[i]=np.mean(np.array([D[j] for j in np.arange(i-3,i+4)]))
    else:
        D_N[i]=D[i]

DD = np.array([(D_N[i+1]-D_N[i-1])/(Z[i+1]-Z[i-1]) for i in np.arange(1,Z.shape[0]-1)])
DD = np.concatenate((np.array([0.]),DD))
DD = np.concatenate((DD,np.array([0.])))
f = interp1d(Z,D_N)
fp = interp1d(Z,DD)
if __name__ == "__main__":
    plt.plot(Z,D_N,color='r')
    plt.plot(Z,D,color='b')
    plt.scatter(Z,D)
    plt.scatter(Z,D_N)
    plt.plot(Z,DD,color='green')
    plt.scatter(Z,DD,color='green')
    plt.plot(Z,f(Z),color='k')
    plt.plot(Z,fp(Z),color='k')
    plt.show()
