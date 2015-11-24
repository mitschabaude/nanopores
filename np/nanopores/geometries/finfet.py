from ..tools.box import Box, union

# geo parameters
lb = 7.
wb = 7.
lw = 80. #length of the fin
hw = hb = 12.
ww = 5.
t1 = 1.5 #thickness of the oxide layer
t2 = 3. #thickness of the gate
ll = 10. #the length of the layers
rdop = 1.
Ndop = 2

# build geometry
source = Box(center=[-(lw + lb)/2, 0, 0], l=lb, w=wb, h=hb)
drain = Box(center=[(lw + lb)/2, 0, 0], l=lb, w=wb, h=hb)
fin = Box(center=[0, 0, 0], l=lw, w=ww, h=hw)
oxide = Box(center=[0, 0, t1/2], l=ll, w= ww + 2*t1, h= hw + t1)
gate = Box(center=[0, 0, (t1 + t2)/2], l=ll, w= ww + 2*(t1 + t2), h= hw + (t1 + t2))

finfet = fin | source | drain | oxide | gate

finfet.addsubdomains(
    sourcendrain = source | drain,
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


