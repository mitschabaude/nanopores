# (c) 2017 Gregor Mitscha-Baude
from simulation import *

dt = .0002
L = 1.
rion = .15
rwater = .15
Nw = 30

cation = Particle([.5*L, 0., .5*L], a=rion, charge=1, color="red")
anion = Particle([-.5*L, 0., -.5*L], a=rion, charge=-1, color="blue")

cationp = partial(random_particle, L, a=rion, charge=1, color="red")
anionm = partial(random_particle, L, a=rion, charge=-1, color="blue")
water = partial(random_particle, L, a=rwater, charge=0, color="#cccccc")

#pcoll = ParticleCollection([cation, anion])
pcoll = ParticleCollection([])
pcoll.add(cationp, 2)
pcoll.add(anionm, 2)
#pcoll.add(water, Nw)

ax, patches, boundaries = setup_box(L)
ani = panimate(ax, dt, pcoll, patches, boundaries, frames=1800)
maybe_save(True, ani, "two_ions_only.mp4")
