from oct2py import Oct2Py
PATH = "collocation_methods/"

# general stuff (wrapper of matlab function)
def collocation(N, order, a0, l, w, h, r):
    octave = Oct2Py()
    octave.addpath(PATH)
    A, W = octave.Collocation_method_sparse_grids(a0, order, l, h, w, r, N)
    octave.exit()
    return A, W

# test
if __name__ == "__main__":
    N = 16
    A, W = collocation(N=N, order=2, a0=[0.,0.,0.], l=8., w=2., h=2., r=0.5)
    print(A.shape)
    print(W.shape)
    print("First sample dopant positions: (N=%d)" %N)
    for i in range(N):
        print(" #", i+1, " ", A[3*i:3*i+3, 0])
    print("Weights:", W)
        

# finfet stuff
import numpy
from nanopores.geometries.finfet import lb, wb, hb, lw, rdop
epsilon = 0.0
rdop1 = rdop + epsilon

def dopants(Ndop, order=2):
    aleft = [-lb - lw/2., -wb/2., -hb/2.]
    aright = [lw/2., -wb/2., -hb/2.]
    
    Aleft, weights = collocation(N=Ndop, order=order, a0=aleft, l=lb, w=wb, h=hb, r=rdop1)
    Aright, weights = collocation(N=Ndop, order=order, a0=aright, l=lb, w=wb, h=hb, r=rdop1)
    
    A = numpy.concatenate([Aleft, Aright])
    Ncol = A.shape[1]
    Ndop = A.shape[0]/3
    dops = []
    
    for i, a in enumerate(A.T):
        dops.append([None]*Ndop)
        for j in range(Ndop):
            dops[i][j] = list(a[3*j:3*j+3])
            
    if isinstance(weights, float):
        weights = [[weights]]
    return dops, list(weights[0])
    
# test
if __name__ == "__main__":
    from nanopores import add_params
    add_params(
        N = 4,
        order = 2,
    )
    dops, weights = dopants(N, order)
    
    for i, sample in enumerate(dops):
        print("Sample # %d:" %i)
        for j, dop in enumerate(sample):
            print("  dopant #%d: %s" %(j, dop))
    
    print("Weights:")
    print(weights)
