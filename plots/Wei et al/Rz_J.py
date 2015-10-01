import matplotlib.pyplot as plt
import numpy as np

# The parameters are according to Fig 1f

slabel = ["N",
          "hmin",
          "distance of electrodes [nm]",
          "current [nA]",
          ]

N_hmin_Rz_Javg = [[47976, 1.806876739528059e-11, 1.5000000000000002e-07, 28532.397863148115], [46700, 1.7921705736806668e-11, 3e-07, 27039.155963142577], [45928, 1.7784629594654474e-11, 4.5000000000000003e-07, 25727.085146792357], [46868, 1.7784629594654474e-11, 6e-07, 24528.456442395567], [48692, 1.7784629594654474e-11, 7.5e-07, 23438.014264718397], [49968, 1.7784629594654474e-11, 9.000000000000001e-07, 22499.9405240072]]

qoilist = [
    N_hmin_Rz_Javg,
]

filestr = ["2015-06-17_13h22.json"]

sinfo = {"comptime": 604.1777172088623, 
         "datetime": "2015-06-17 13:22:02.526393", 
         "geo_name": "W_2D_geo", 
         "params": {
             "bV": 0.2, 
             "bulkcon": 1000.0, 
             "lsam": 3.0000000000000004e-09, 
             "r0": 1e-08, 
             "x0": None, }
     }

sleg = [
    "pnps",
]

for i in range(len(qoilist)):
    qoi = qoilist[i]
    N=[]
    hmin=[]
    Rz=[]
    J=[]
    
    for j in range(len(qoi)):
        N.append(qoi[j][0])
        hmin.append(qoi[j][1])
        Rz.append(qoi[j][2]*2*1e9)
        J.append(qoi[j][3])        

    # #ran = slice(0,len(hmin),1)
    plt.figure(1)
    plt.plot(Rz, J, 'x-', label=sleg[i] )
    plt.xlabel(slabel[2])
    plt.ylabel(slabel[3])
    legend = plt.legend() #loc='lower right')
    
    # plt.figure(2)
    # plt.plot(N, J, 's-', label=sleg[i] )
    # plt.xlabel(slabel[3])
    # plt.ylabel(slabel[2])
    # legend = plt.legend(loc='upper left')
plt.show()
