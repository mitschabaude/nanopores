# rx = 15 probably much too low
ry = [700, 600, 500, 400, 300, 200, 150, 100, 75, 50, 25, 15, 11]
bV = [0.0397, 0.0435, 0.0480, 0.0536, 0.0588, 0.0698, 0.0756, 0.0813, 0.0853, 0.0897, 0.0946, 0.0967, 0.09837]
i = [194, 213, 235, 262, 288, 341, 368, 451, 472, 497, 524, 534, 539]

# rx = ry this should be enough but scales badly..
ry = [100, 80, 60, 40, 30, 20, 15]
bV = [0.0957, 0.0957, 0.0958, 0.0960, 0.09614, 0.09648, 0.09700]
i = [466, 467, 466, 464, 467, 468, 468, 470]

# at this point i have no claim at all that the membrane potential is effectively reduced by a factor 10 or more :(
# but we can still reduce the diffusion constant by a factor 10 :)

import matplotlib.pyplot as plt

fname = "effbV.eps"

#plt.loglog(ry,bV,'s-')
#plt.semilogx(ry,bV,'s-')
plt.plot(ry,bV,'s-')
plt.xlabel("Electrode distance [nm]")
plt.ylabel("Eff. transmembrane voltage [mV]")
plt.show()

