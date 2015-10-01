#!/usr/bin/env python
''' create mesh and midpoints file for force calculations '''

from nanopores import *
from dolfin import *
from numpy import *

from protocol import Data
import os, sys

DIR = DATADIR + "/sim"
simname = "sim"

if not os.path.exists(DIR):
    os.makedirs(DIR)
    
if len(sys.argv) > 1:
    simname = sys.argv[1]

FILENAME = DIR + "/" + simname + ".dat"
ndigits = 9+3

def create_data_file(points):
    points = array(points)
    N = len(points)

    data = Data(FILENAME, N=N, overwrite=True)
    data["r"] = zeros(N)
    data["z"] = points

    return data.filename

def drange(start, stop, step):
    r = start
    while min(start,stop) <= r <= max(start,stop):
        yield r
        r += step

nm = 1e-0
step = 0.2*nm
Rz = 12*nm

z = list(drange(-Rz, Rz, step))

print "Created data file in", create_data_file(z)
print "Number of data points:", len(z)
