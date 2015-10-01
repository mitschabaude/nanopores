import matplotlib.pyplot as plt
import numpy as np


sollabel = ["SiN Surface Charge [mC/m^2]",
            "Effective driving force [pN]"]
sol32k = [np.arange(0,-70,-10),
          np.array([11.1003771811, 10.2954075218, 9.36374745865, 8.9678418328,
                 7.7674555794, 7.58699033813, 6.49438780848])]
sol64k = [np.arange(0,-70,-10),
          np.array([10.8803657169, 10.0334827643, 9.05291908243, 8.32836358782,
                    7.25693958528, 6.60970424053, 5.78215811963])]
pre = [np.arange(0,-70,-10),
       np.array([12.2, 11.15, 10.15, 9.1, 8.15, 7.2, 6.25])]

fig, subp = plt.subplots()
subp.plot(sol64k[0], sol64k[1], 'rx-', label="adaptive FEM solution with 64k cells",)
#subp.plot(sol32k[0], sol32k[1], 'bo-',)
subp.plot(pre[0], pre[1], 'ks-', label="Effective driving force on DNA... (PRE)",)

subp.axis([-60, 0, 5, 13])
plt.xlabel(sollabel[0])
plt.ylabel(sollabel[1])
#plt.title("%s over %s" %(sollabel[0], sollabel[1]))
legend = subp.legend(loc='upper center')

plt.show()


                 
