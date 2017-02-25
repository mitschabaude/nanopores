#import matplotlib
#matplotlib.use("Agg")
#import matplotlib.pyplot as plt
from random import gauss
from math import sqrt
import numpy as np
import sys

taun = sys.argv[1]
tau = float(sys.argv[1])
runs = int(sys.argv[2])
eps = float(sys.argv[3])
coeff = sqrt(tau)

data = np.load('data2'+taun+'.npy')

print 'tau = %.8f'%tau
for i in range(runs):
    y_old = 0.
    time = 0.
    maxtime = 5.
    attempt = 0
    out=0.
    while time<=maxtime:
        xi = gauss(0.,1.)
        y_new = y_old + coeff*xi
        if abs(y_new)>eps:
            if np.sign(out)!=np.sign(y_new):
                attempt+=1
            out=y_new
        time+=tau
        y_old = y_new
    data = np.append(data,attempt)
    if i%1000 == 0:
        print '%i out of %i'%(i,runs)
np.save('data2'+taun,data)
