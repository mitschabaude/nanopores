from ..tools.box import Box, union

# geo parameters
lb = hb = 6.
wb = 8.
lw = 30.
hw = ww = 4.
t1 = .5
t2 = 2.
ll = 2.
rdop = 1.
Ndop = 2

# build geometry
leftblock = Box(center=[-(lw + hb)/2, 0, (hb-hw)/2], l=lb, w=wb, h=hb)
rightblock = Box(center=[(lw + hb)/2, 0, (hb-hw)/2], l=lb, w=wb, h=hb)
wire = Box(center=[0, 0, 0], l=lw, w=ww, h=hw)
layer1 = Box(center=[0, 0, t1/2], l=ll, w= ww + 2*t1, h= hw + t1)
layer2 = Box(center=[0, 0, (t1 + t2)/2], l=ll, w= ww + 2*(t1 + t2), h= hw + (t1 + t2))

finfet = wire | leftblock | rightblock | wire | layer1 | layer2

finfet.addsubdomains(
    blocks = leftblock | rightblock,
    wire = wire,
    layer1 = layer1 - wire,
    layer2 = layer2 - layer1,
)

finfet.addboundaries(
    gate = layer2.boundary("top", "front", "back"),
    leftb = leftblock.boundary("left"),
    rightb = rightblock.boundary("right"),
)



