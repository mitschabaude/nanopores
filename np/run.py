from random import gauss
from math import sqrt, pi
import math
import numpy as np
from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
import sys
from calculateforce import *
from aHem_array import *
F = calculateforce(clscale=10., tol=5e-1) # 6. 5e-3
# def F(vec):
#     return [0,0,-0.01]
def radius(x,y):
    return sqrt(x**2+y**2)

def normal(ax,ay,bx,by,px,py):
	AP2=(ax-px)**2+(ay-py)**2
	BP2=(bx-px)**2+(by-py)**2
	AB2=(ax-bx)**2+(ay-by)**2
	AB=sqrt(AB2)
	c = (AP2-BP2+AB2)/(2*AB)
	if c>0. and c<AB:
		if AP2<=c**2:
			return 0.
		else:
			return sqrt(AP2-c**2)
	else:
		return 100.

def distance_to_surface(rad,z):
	if z>3.:
		return 100.
	elif rad>8.:
		return 100.
	else:
		size=X_aHem.shape[0]
		D = np.zeros(size)
		E = np.zeros(size)
		for index in range(size):
			D[index] = radius(rad-X_aHem[index][0],z-X_aHem[index][2])
			E[index] = normal(X_aHem[(index-1)][0],X_aHem[(index-1)][2],X_aHem[(index)][0],X_aHem[(index)][2],rad,z)
		return min([min(D),min(E)])

def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])

geo = geo_from_xml("aHem")
indicator_ahem_geo = geo.indicator("ahem",callable=True)
indicator_molecule = geo.indicator("molecule",callable=True)
indicator_porecenter_geo = geo.indicator("porecenter",callable=True)
indicator_membrane_geo = geo.indicator("membrane",callable=True)
indicator_poretop_geo = geo.indicator("poretop",callable=True)
def indicator_membrane(vec): #"infinite" large membrane
    x, y, z = vec[0], vec[1], vec[2]
    if z>60.:
    	return 0
    elif radius(x,y)>=30.:
        return indicator_membrane_geo(argument(10,0,z))
    else:
        return indicator_membrane_geo(vec)
def indicator_ahem(vec):
	x, y, z = vec[0], vec[1], vec[2]
	if radius(x,y)>=30. or z>60.:
		return 0
	else:
		return indicator_ahem_geo(vec)
def indicator_porecenter(vec):
    x, y, z = vec[0], vec[1], vec[2]
    if radius(x,y)>50.:
    	return 0
    elif z>0.:
    	return 0
    else:
    	return indicator_porecenter_geo(vec)
def indicator_poretop(vec):
    x, y, z = vec[0], vec[1], vec[2]
    if radius(x,y)>50.:
    	return 0
    elif z>0.:
    	return 0
    else:
    	return indicator_poretop_geo(vec)

def oor(x,y,z):
    return radius(x,y)>500 or z>500


kb=1.3806488e-23 #boltzmann [J/K]
T= 293 #temp [K]
visc=1e-3 #[Pa*s]
D=(kb*T)/(6*pi*0.5*visc)*1e9  #diffusion[m^2/s]
gamma=(6*pi*0.5*visc) #friction [kg/s]
tau=0.1 #-1 # [ns]
steps=5e6 # should be 5e7
C=1/(gamma)*tau
coeff=sqrt(2*D*1e9*tau)

counter = np.array([0,0,0,0,0]) # ahem, molecule, poretop, membrane, bulk
EXIT_X, EXIT_Y, EXIT_Z = np.array([]), np.array([]), np.array([])
Range =range(1)# range(300)
for index in Range:
    print str(index)+" out of "+str(len(Range))
    X=np.zeros((steps))
    Y=np.zeros((steps))
    Z=np.zeros((steps))

    Z[0] = 2.
    i=0
    while i<=steps-2:
        a=i/(steps-2)*100
        rad=radius(X[i],Y[i])
        print
        print ">>>>>>>>>>>>>>>> progress %0.1f percent" %a
        print 'RAD = %.3f , Z = %.3f '%(rad, Z[i])
        sys.stdout.write("\033[F") # Cursor up one line
        sys.stdout.write("\033[F") # Cursor up one line
        sys.stdout.write("\033[F") # Cursor up one line
        xi_x=gauss(0,1)
        xi_y=gauss(0,1)
        xi_z=gauss(0,1)
        if indicator_poretop(argument(X[i],Y[i],Z[i]))==1: #Targetmolecule in Pore ==> diffusion damp
            xi_x *= 0.316
            xi_y *= 0.316
            xi_z *= 0.316
        X[i+1]=X[i] + coeff*xi_x + C*F(argument(X[i],Y[i],Z[i]))[0]
        Y[i+1]=Y[i] + coeff*xi_y + C*F(argument(X[i],Y[i],Z[i]))[1]
        Z[i+1]=Z[i] + coeff*xi_z + C*F(argument(X[i],Y[i],Z[i]))[2]
        if indicator_membrane(argument(X[i+1],Y[i+1],Z[i+1]))==1:
            exit_x = X[i+1]
            exit_y = Y[i+1]
            exit_z = Z[i+1]
            counter[3] += 1
            i+=1
            break
        elif indicator_ahem(argument(X[i+1],Y[i+1],Z[i+1]))==1:
            exit_x = X[i+1]
            exit_y = Y[i+1]
            exit_z = Z[i+1]
            counter[0] += 1
            i+=1
            break
        elif indicator_porecenter(argument(X[i+1],Y[i+1],Z[i+1]))==1:
            exit_x = X[i+1]
            exit_y = Y[i+1]
            exit_z = Z[i+1]
            counter[2] += 1
            i+=1
            break
        # elif indicator_molecule(argument(X[i+1],Y[i+1],Z[i+1]))==1:
        #     exit_x = X[i+1]
        #     exit_y = Y[i+1]
        #     exit_z = Z[i+1]
        #     counter[1] += 1
        #     i+=1
        #     break
        if oor(X[i+1],Y[i+1],Z[i+1]):
            i=steps-2
        i+=1
    if i>steps-2:
        counter[4] += 1
    else:
        EXIT_X = np.append(EXIT_X, np.array([exit_x]))
        EXIT_Y = np.append(EXIT_Y, np.array([exit_y]))
        EXIT_Z = np.append(EXIT_Z, np.array([exit_z]))
np.save('exit_x',EXIT_X)
np.save('exit_y',EXIT_Y)
np.save('exit_z',EXIT_Z)
np.save('counter',counter)

import plot