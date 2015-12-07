import nanopores
import numpy
import dolfin
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

dolfin.tic()

bmV = numpy.linspace(-99.9, 99.9, 8)
I = []

# IV ohmic and diodic fit taken from:
# Rectification of the Current in alpha-Hemolysin Pore Depends on the Cation Type:
# The Alkali Series Probed by Molecular Dynamics Simulations and Experiments
# Here, the currents were calculated by using the following parameters
# for equation I = bV/(r+R) + r*I0/(r+R)*(numpy.exp(bV-R*I/V0)-1)

sigmab = 108  # [mS/cm]
I0 = 24.  # [pA]
V0 = 4.  # [mV]
Gminf = 0.75  # 1/(R+r) [nS]
Gpinf = 1.15  # 1/R [nS]
# Gpinf/sigmab = 1.0 Angstrom

R = 1/Gpinf  # [G Ohm]
r = 1/Gminf - R  # [G Ohm]

for amV in bmV:
    func = lambda i: i - amV/(r+R) - r*I0/(r+R)*(numpy.exp((amV-R*i)/V0)-1)
    I_guess = amV/(R+r)*2  # [V/G Ohm *1e-3 = pA]
    I_sol = fsolve(func, I_guess)[0]
    I.append(I_sol)  # [pA]
    # print I_guess, I_sol

print "\nbmV = ", bmV
print "I = ", I

### SIMULATION

# or use existing .dat file with uid
uid = None

# Simulation Args
bV = bmV*1e-3
args = dict(
    nproc = 1,
    plot = "bV",
    outputs = ["Javgbtm", "Javgtop"],

    domscale = 1,
    #r0 = 5.,
    #z0 = 10.,

    bV = list(bV),
    bulkcon = 1000.,
    ahemqs = list(numpy.linspace(0.0, -0.08, 3)),
    rDPore = 0.22, # list(numpy.linspace(0.2, 0.3, 5)),
    #rTarget = [3e-10, 10e-10, 1e-9]

    clscale = 6.,
    skip_stokes = True, #[True, False],
    iterative = True,)

if not uid:
    sim = nanopores.simulate("ahemIV",**args)
else:
    from nanopores.tools.protocol import *
    from nanopores import DATADIR
    savedir = DATADIR + "/sim/stamps/"
    data = Data(savedir+"result"+uid+".dat")
    datadict = data.data
    # print datadict.keys()

    Isim = datadict["Javgtop"]
    # specific, assume data files corresponds to this simulation inputs
    bVsim = args["bV"]
    qssim = numpy.asarray(args["ahemqs"])*1e3
    # print "Isim: ", Isim
    # print "bVsim: ", bVsim

    # regroup
    Isim0 = Isim
    nx = len(bVsim)
    ny = range(int(float(len(Isim0))/nx))
    Isim0p = [Isim0[slice(i*nx,(i+1)*nx)] for i in ny]

    plt.figure()
    for i in ny:
        plt.plot(bV, Isim0p[i], '-x', label='PNPS '+str(qssim[i])+' mC/m$^2$')

# # get experimental data from csv file
# # Experimental Values taken from:
# # Properties of Bacillus cereus hemolysin II: A heptameric transmembrane pore
# csvfile = 'data/ahemIV/ahemIV'
# IV = numpy.genfromtxt(csvfile+'.csv', delimiter=',')
# Vexp = -IV[1:-1,0]*1e-3
# Iexp = -IV[1:-1,1]
# plt.plot(Vexp, Iexp, '-ks', label='experimental')

plt.plot(bV, I, '-^', label='ICR (+ diode)')
plt.plot(bV, bmV/(R+r), '-v', label='ohmic')
plt.xlabel("V [V]")
plt.ylabel("I [pA]")
title = ("relative diffusivity in pore: "+str(args["rDPore"]) if not hasattr(args["rDPore"], '__iter__') else "ahem surface charge: "+str(args["ahemqs"]))
plt.title(title)
plt.legend(loc=0,)
plt.grid('on')
plt.show()

print "\nTotal Time: ", dolfin.toc()
