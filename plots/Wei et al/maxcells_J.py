import matplotlib.pyplot as plt
import numpy as np

# This plot reproduces result Fig S2b in Wei Supplementary info

slabel = ["N",
          "hmin",
          "maxcells",
          "current [nA]",
          ]

N_hmin_maxcells_Javg_pnp = [[6082, 2.3555347807025036e-10, 5000.0, 26836.108000983673], [12073, 1.1941179949055674e-10, 14000.0, 28570.31513715205], [20187, 5.2602573906724294e-11, 23000.0, 28560.62026720029], [30275, 3.113809382745503e-11, 32000.0, 28550.463056925648], [30275, 3.113809382745503e-11, 41000.0, 28550.444312774547], [46743, 1.7485622762228622e-11, 50000.0, 28531.7561122048]]
#filestr = '2015-06-17_15h03.json'

N_hmin_maxcells_Javg_pnps = [[6082, 2.3555347807025036e-10, 5000.0, 26830.84952099491], [12079, 1.1941179949055674e-10, 14000.0, 28562.511042714155], [20193, 5.2602573906724294e-11, 23000.0, 28552.83734644027], [30305, 3.113809382745503e-11, 32000.0, 28542.659642088085], [30305, 3.113809382745503e-11, 41000.0, 28542.64095924403], [46777, 1.7485622762228622e-11, 50000.0, 28523.636139677837]]
#filestr = '2015-06-17_15h01.json'

qoilist = [
    N_hmin_maxcells_Javg_pnp,
    N_hmin_maxcells_Javg_pnps,
]

sleg = [
    "pnp - sam 1f",
    "pnps - sam 1f",
]

for i in range(len(qoilist)):
    qoi = qoilist[i]
    N=[]
    hmin=[]
    vparam=[]
    J=[]
    
    for j in range(len(qoi)):
        N.append(qoi[j][0])
        hmin.append(qoi[j][1])
        vparam.append(qoi[j][2])
        J.append(qoi[j][3])        

    # #ran = slice(0,len(hmin),1)
    plt.figure(1)
    plt.plot(N, J, 'x-', label=sleg[i] )
    plt.xlabel(slabel[0])
    plt.ylabel(slabel[3])
    legend = plt.legend(loc='lower right')
    
    plt.figure(2)
    plt.semilogx(hmin, J, 's-', label=sleg[i] )
    plt.xlabel(slabel[1])
    plt.ylabel(slabel[3])
    legend = plt.legend(loc='upper right')
plt.show()
