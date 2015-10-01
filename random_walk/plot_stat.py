import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
infile = open('file.txt', 'r')
data = infile.read()
infile.close()
data = data.split()
i=1
while i<=len(data):
    data[i-1]=float(data[i-1])
    i+=1
Y=np.zeros(len(data))
X=np.array(data)
plt.scatter(X,Y)
plt.show()
