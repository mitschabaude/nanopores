from scipy import ndimage
import nanopores.geometries.pughpore as pughpore
import nanopores
import numpy as np
import matplotlib.pyplot as plt

up = nanopores.user_params(pughpore.params, k=3)

R = up.R
H = up.H

def z_func(x, y):
    return (1 - (x ** 2 + y ** 3)) * np.exp(-(x ** 2 + y ** 2) / 2)

X, Y, Z = np.load('X.npy'), np.load('Y.npy'), np.load('ZJ4.npy')
Z=ndimage.rotate(Z,0)

im = plt.imshow(Z, cmap=plt.cm.viridis, extent=(-R/2., R/2., H/2., -H/2.))  
plt.colorbar(im)  

plt.tight_layout()
plt.show()
