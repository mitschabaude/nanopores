from ..tools.box import Box, union
from random import random

# geo parameters
lb = 7.
wb = 7.
lw = 80. #length of the fin
hw = hb = 7.
ww = 5.
t1 = 1.5 #thickness of the oxide layer
t2 = 3. #thickness of the gate
ll = 10. #the length of the layers
rdop = 1.
Ndop = 12

# build geometry
source = Box(center=[-(lw + lb)/2, 0, 0], l=lb, w=wb, h=hb)
drain = Box(center=[(lw + lb)/2, 0, 0], l=lb, w=wb, h=hb)
fin = Box(center=[0, 0, 0], l=lw, w=ww, h=hw)
oxide = Box(center=[0, 0, t1/2], l=ll, w= ww + 2*t1, h= hw + t1)
gate = Box(center=[0, 0, (t1 + t2)/2], l=ll, w= ww + 2*(t1 + t2), h= hw + (t1 + t2))

finfet = fin | source | drain | oxide | gate

finfet.addsubdomains(
    source = source,
    drain = drain,
    fin = fin,
    oxide = oxide - fin,
    gate = gate - oxide,
)

finfet.addboundaries(
    gateb = gate.boundary("top", "front", "back"),
    sourceb = source.boundary("left"),
    drainb = drain.boundary("right"),
    crossl = fin.boundary("left"),
    crossr = fin.boundary("right"),
)

finfet.synonymes = dict(
    sourcendrain = {"source", "drain"}
)

# add parameters (this should include params needed by physics module)
finfet.params = dict(
    rdop = rdop,
    lscale = 1e9,
    length = lw + 2.*lb,
    lfin = lw
)


# define dopant creation function here because this needs the box sizes
def dopants(Ndop=Ndop):
    # create two chunks of random vars in [0, 1]**3
    X = [[random() for i in range(3)] for i in range(Ndop)]
    X1 = [[2*x[0] - 1., x[1], x[2]] for x in X if x[0] > .5]
    X2 = [[2*x[0], x[1], x[2]] for x in X if x[0] <= .5]
    # affine transform to source and drain
    def affine(X, box, R):
        return [[ai + R + t*(bi-ai-2.*R) for t, ai, bi in zip(x, box.a, box.b)] for x in X]
    X1 = affine(X1, source, rdop)
    X2 = affine(X2, drain, rdop)
    return X1 + X2

