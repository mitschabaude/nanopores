# enhanced Stokes drag on spherical particle in cylindrical tube

# r ... particle radius
# R ... cylinder radius

# drag ... reduced stokes drag

# dimensionless parameters:
# lam = r/R in [0,1]
# mu  = r/(R-r) in [0,inf]
# rho = (R-r)/r = 1/lam - 1 in [0,inf]

# empirical law:
# drag = exp(2*r/(R-r)) = exp(2*lam/(1-lam)) = exp(2*mu) = exp(2/rho)

# get experimental data from csv file
from numpy import genfromtxt, exp
TT = genfromtxt('stokes_drag_table_rdivR.csv', delimiter=',')
lam = TT[:,0] # lam = r/R
drag0 = TT[:,1]
drag = [exp(2.*x/(1.-x)) for x in lam]

D0 = [1./x for x in drag0]
D = [1./x for x in drag]

rho = [1./x - 1. for x in lam[::-1]]
drag0rho = drag0[::-1]
dragrho = [exp(2./x) for x in rho]

import matplotlib.pyplot as plt
'''
#fname = "IV_c0.1_Ry1000nm.eps"
plt.plot(lam,drag0, 's-', label="Paine, Scherr")
plt.plot(lam,drag, 'x-', label="Mitscha-Baude")
plt.xlabel("r/R")
plt.ylabel("drag/drag0")
plt.legend(loc='lower right')

#plt.savefig(fname, bbox_inches='tight')
plt.show()
'''
#fname = "IV_c0.1_Ry1000nm.eps"
plt.plot(lam,D0, 's-', label="Paine, Scherr")
plt.plot(lam,D, 'x-', label="Mitscha-Baude")
plt.xlabel("r/R")
plt.ylabel("D/D0")
plt.legend(loc='upper right')
'''
#fname = "IV_c0.1_Ry1000nm.eps"
plt.loglog(rho,drag0rho, 's-', label="Paine, Scherr")
plt.loglog(rho,dragrho, 'x-', label="Mitscha-Baude")
plt.xlabel("(R-r)/r")
plt.ylabel("drag/drag0")
plt.legend(loc='upper right')
'''
#plt.savefig(fname, bbox_inches='tight')
plt.show()

