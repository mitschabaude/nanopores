import matplotlib.pyplot as plt
import numpy as np

slabel = ["N", "hmin", "Surface Charge [C/m^2]", "current [pA]"]

filestr = "2015-06-17_12h58.json"

N_hmin_Membraneqs_Javg = [[43344, 4.3235002624703605e-12, 0.0, 28160.214315528356], [50289, 1.748562276226755e-11, -0.02, 28444.288825520227], [47114, 1.8786630763613166e-11, -0.04, 28660.70559658548], [46776, 1.7784629594654474e-11, -0.06, 28962.22167173465], [45168, 1.8786630763613166e-11, -0.08, 29279.170169005425], [44278, 1.961385267047395e-11, -0.1, 29608.379007002073], [43233, 1.961385267047395e-11, -0.12, 29963.566932733414], [42733, 1.961385267047395e-11, -0.14, 30325.170914428432], [43001, 1.9613852971815267e-11, -0.16, 30667.11857003202], [42934, 1.9613852670527288e-11, -0.18, 31004.312936641367], [41546, 1.9613852670509333e-11, -0.2, 31365.22793016157]]

qoilist = [N_hmin_Membraneqs_Javg, ]

sleg = [
    "pnps - sam - 50k cells",
]

for i in range(len(qoilist)):
    qoi = qoilist[i]
    N=[]
    hmin=[]
    qs=[]
    J=[]
    
    for j in range(len(qoi)):
        N.append(qoi[j][0])
        hmin.append(qoi[j][1])
        qs.append(qoi[j][2])
        J.append(qoi[j][3])

    # #ran = slice(0,len(hmin),1)
    plt.figure(1)
    plt.plot(qs, J, 'x-', label=sleg[i] )
    plt.xlabel(slabel[2])
    plt.ylabel(slabel[3])
    legend = plt.legend(loc='lower right')
    
    # plt.figure(2)
    # plt.plot(N, J, 's-', label=sleg[i] )
    # plt.xlabel(slabel[3])
    # plt.ylabel(slabel[2])
    # legend = plt.legend(loc='upper left')

plt.show()        
