#!/usr/bin/env python

''' calculations with N points (passed by user) '''

from nanopores import *
from dolfin import *
from numpy import *
from protocol import Data, get_subindices, get_remaining
import sys, importlib, subprocess
from argparse import ArgumentParser
from calculate_forces import calculate_forces, calculate_forces2D

# which forces you want to write to data file: (subset of pnps.functionals)
#fkeys = {"Javgtop","Javgctr","Javgbtm","Fp[0]","Fp[1]","Fp[2]","Fshear[0]",
#        "Fshear[1]","Fshear[2]","Fbare[0]","Fbare[1]","Fbare[2]"}
fkeys = {"Javgtop","Javgctr","Javgbtm","Fp0","Fp1","Fp2","Fshear0",
        "Fshear1","Fshear2","Fbare0","Fbare1","Fbare2","Fbarevol0","Fbarevol1","Fbarevol2"}

def prettytime(t):
    return ":".join(["%02d"%int(k) for k in (t//3600, (t%3600)//60, (t%3600)%60)])
    
def run_with_output(args, verbose=True):
    ''' run subprocess, display output and return last line of output '''
    # start process
    process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # continually check output until process terminates
    output = save = process.stdout.readline()
    while ((output != '') or (process.poll() is None)):
        if output:
            if verbose:
                print "    ",output,
            save = output
        output = process.stdout.readline()
    # return exit code and last line of output
    return process.poll(), save.strip()

t0 = Timer("Total time")
DIR = DATADIR + "/sim"
FILENAME = DIR + "/sim.dat"

# default values for user arguments
N = 1
CLSCALE = 4.0

# parse user arguments
parser = ArgumentParser()
parser.add_argument('-N','--nmol', default=N, type=int, help='number of molecule positions')
parser.add_argument('-p','--pid', default="", help='process ID for .out .msh .xml files')
parser.add_argument('-v','--verbose', action="store_true", help='if given, calculation output is printed')
parser.add_argument('-c','--clscale', default=CLSCALE, type=float, help='mesh width scale')
parser.add_argument('-2','--dim2', action="store_true", help='do 2D calculation (default is 3D)')
args, unknown = parser.parse_known_args()
print vars(args)
N = args.nmol
pid = args.pid
CLSCALE = args.clscale
dim = 2 if args.dim2 else 3
if args.dim2:
    calculate_forces = calculate_forces2D

# We assume that FILENAME has been created by ./create_mesh.py
t1 = Timer("Data access time")
data = Data(FILENAME)
subdata, uid = data.get_points(N, read=False)
N = subdata["i"].size # actual number of returned points
rem = get_remaining(data.data)

print "\nAccessed file %s, got %s vertices." %(FILENAME, str(N))
print "Data access time:",t1.stop()
print "%s vertices remaining." %rem

# from here on, data is not accessed until the very end
# this means, other processes can safely communicate with data file now

# initialize/overwrite forces and get coordinates
for key in fkeys:
    subdata[key] = zeros(N)

index = subdata["i"]
points = zip(subdata["r"], subdata["z"])
calculated = []


# GOGOGO
for i in range(N):
    r,z = points[i]
    x0 = [r, 0.0, z]
    
    #if i>0: print "\x1b[A"*6 + "\r",
    print "\nMolecule loop",i+1,"of",N
    print "x0 =",x0
    print "Calculating..."
    try:
        t = Timer('Force calculation')
        
        # start calculation for one molecule
        forces = calculate_forces(x0, pid=pid, clscale=CLSCALE,
                                  refinement=True)
        #excode, out = run_with_output(["python", "calculate_forces.py",
        #     repr(x0), pid, str(CLSCALE), str(dim)], verbose=args.verbose)
        
        # read forces from last line of output and save
        #forces = eval(out)
        
        for key in forces:
            subdata[key][i] = forces[key]    
        calculated.append(i)
        
        # print timing / remaining time
        if i < N-1:
            #print "\x1b[A\r",
            s = t.stop()
            print "\nLoop time:",prettytime(s)
            print "Projected remaining time:",prettytime(s*(N-1-i))
    #except (KeyError, ) as e:        
    except (Exception, KeyboardInterrupt) as e:
        #print "\x1b[A\r",
        print "\nSome error occured, no results were saved."
        print ": ".join([type(e).__name__, str(e)])

s = t0.stop()    
print "\nTotal time for %s vertices:" %N, prettytime(s)        

subdata = get_subindices(subdata, calculated)

# write back to data
t1.start()
remaining = data.write_values(subdata, uid)
print "Data access time:",t1.stop()

print "\nWrote results for %s vertice(s) to data file.\n%s vertice(s) remaining." %(len(calculated), remaining)
print

