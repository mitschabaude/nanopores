# (c) 2017 Gregor Mitscha-Baude
import nanopores
from nanopores.models.randomwalk import (get_pore, RandomWalk, video, Ball,
                                         plt, histogram)
import matplotlib as mpl
mpl.verbose.set_level("helpful")
params = nanopores.user_params(
    # general params
    geoname = "alphahem",
    dim = 2,
    rMolecule = .25,
    h = 1.,
    Nmax = 4e4,
    Qmol = 1.,
    bV = -0.4,
    # random walk params
    N = 100, # number of (simultaneous) random walks
    dt = 10., # time step [ns]
    walldist = 2., # in multiples of radius, should be >= 1
    margtop = 2.,
    margbot = .5,
    rstart = 0.1
)

pore = get_pore(**params)
rw = RandomWalk(pore, **params)
#receptor = Ball([1., 0., -30.], 8.)
#rw.add_domain(receptor, exclusion=True, walldist=1.,
#              binding=True, eps=1., t=1e6, p=0.1)
#rw.add_wall_binding(t=1e4, p=0.1, eps=0.1)
videop = nanopores.user_params(save=False, cyl=False, video=False)

if videop.video:
    ani = video(rw, cyl=videop.cyl, save_count=1000)
    if videop.save:
        ani.save(nanopores.HOME + "/presentations/nanopores/ahem.mp4",
                 fps=30, dpi=200,
                 writer="ffmpeg_file",
                 #savefig_kwargs={"bbox_inches":0, "pad_inches":0},
                 extra_args=['-vcodec', 'libx264'])
    else:
        plt.show()
else:
    for t in rw.walk(): pass

if not videop.save:
    histogram(rw, a=-3, b=6)
    plt.show()