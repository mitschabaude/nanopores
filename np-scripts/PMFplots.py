from nanopores import Data, DATADIR, kB, T
from matplotlib.pyplot import *
import os

# integrate vector of function values with uniform trapezoid rule
def antider(F, step):
    U = []
    for f in F:
        if 0 < len(U):
            U.append(U[-1] + 0.5*step*(f+f0))
        else:
            U.append(0.5*step*f)
        f0 = f
    return U

# get git repo
with open("%s/sim/gitrepo"%(DATADIR,), "r") as f:
    gitrepo = f.next().strip()

s = "%s/%s"
folder = s% (gitrepo, "np-sim/simulations")
files = [f for f in os.listdir(folder) if f.startswith("PMF") and f.endswith(".dat")]
#files = ["PMF_r55_DNA0.dat"]
forces = ["Fbarevol2","Fp2","Fshear2"]

nm = 1e-9
step = 0.2*nm

# conversion factor for energy (from pJ)
C = 1e-12/(kB*T) # [kT]
#C = 1e-12*1.44e20 # [kcal/mol]

for file in files:

    # get data from file
    data = Data(s% (folder, file))
    
    # jump to next if calculation is not complete
    if not all(data["status"] == 1): continue
    
    # get force and z coordinates
    F = sum(data[force] for force in forces) 
    Z = data["z"]
    
    # compute PMF in units of [kT]
    # TODO should depend on Z
    U = antider(F, -step*C) #/(kB*T))
    U = map(lambda u : u - U[-1], U)
    #print U[::-1]
    
    figure(1)
    plot(Z, U, '-x', label=file)
    
    figure(2)
    plot(Z, F, '-x', label=file)
    
figure(1)
legend(loc = "lower right")
xlabel("z [nm]")
#ylabel("PMF [kcal/mol]")
ylabel("PMF [kT]")

figure(2)
legend(loc = "upper left")
xlabel("z [nm]")
ylabel("Force [pN]")

show()


