import matplotlib.pyplot as plt
import numpy as np
t=np.load('TAU.npy')*1e-6
a=np.load('A.npy')

plt.scatter(t,a)
ax=plt.gca()
ax.set_xlim([1e-4,500.])
ax.set_ylim([0.,2.5])
ax.set_xscale('log')
ax.invert_yaxis()
ax.set_xlabel('tau_off [ms]')
ax.set_ylabel('A/I_0 [%]')
plt.tight_layout()
plt.show()
