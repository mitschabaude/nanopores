from simulation import *

dt = .001
L = 5
rlarge = 2.1
rsmall = .11
Qlarge = 50
N = 100

Np = (N - Qlarge)/2
Nm = (N + Qlarge)/2

Np = Nm = 66

large = Particle([0.,0.,0.], rlarge, -50., "#9999ff")

cationp = partial(random_particle, L, a=rsmall, charge=1, color="red")
anionm = partial(random_particle, L, a=rsmall, charge=-1, color="blue")

pcoll = ParticleCollection([large])
pcoll.add(cationp, Np)
pcoll.add(anionm, Nm)

ax, patches, boundaries = setup_box(L)
ani = panimate(ax, dt, pcoll, patches, boundaries, frames=1800)
maybe_save(False, ani, "ion_box.mp4")
