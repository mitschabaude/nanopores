# (c) 2016 Gregor Mitscha-Baude
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from nanopores.geometries.pughpore import params as pugh_params
from nanopores import Params

def grid_piecewise1D(nodes, h, N=100, ep=None):
    # compute number of grid points in each section
    # N = 1/scaling * (length1 / h1 + length2 / h2 + ...)
    lengths = np.diff(np.array(nodes))
    h = np.array(h)
    n = lengths/h
    n = np.round(n*N/sum(n))
    # compute each grid
    intervals = zip(nodes[:-1], nodes[1:])
    k = len(lengths)
    grids = []
    # ep = endpoint preference = 0 or 1
    if ep is None:
        ep = [0]*(k-1)
    for i in range(k):
        a, b = intervals[i]
        grid = list(np.linspace(a, b, n[i]+1)[1:-1])
        #print i
        #print grid
        if i == 0 or ep[i-1] == 1:
            grid.insert(0, a)
        if i == k-1 or ep[i] == 0:
            grid.append(b)
        #print grid
        grids.append(grid)
    #print n
    return grids

def lround(x, nd):
    if hasattr(x, "__iter__"):
        return [lround(t, nd) for t in x]
    else:
        return round(x, nd)

def tensor(xy, z, r):
    tensorgrid = []
    for i in range(len(z)):
        # scale xy by radius
        ri = r[i]
        xyz = [(ri*xj, ri*yj, zi) for zi in z[i] for xj, yj in xy]
        tensorgrid.extend(xyz)
    return lround(tensorgrid, 3)


def plot_1Dgrid(z, grids):
    totalgrid = list(set(reduce(lambda a, b: a+b, grids)) - set(z))
    fig = plt.figure("line")
    fig.set_size_inches(8, 1)
    plt.axhline(y=0, color="black", zorder=-10)
    plt.scatter(totalgrid, [0.]*len(totalgrid), color="black")
    plt.scatter(z, [0.]*len(z), color="red")
    plt.xlim(z[0]-1, z[-1]+1)
    plt.axis('off')

def neg(x):
    return [-t for t in x]

def plot_2Dgrid(xy):
    xx = [xi for (xi, yi) in xy]
    yy = [yi for (xi, yi) in xy]
    fig = plt.figure("triangle")
    fig.set_size_inches(4, 4)
    plt.scatter(xx, yy, color="red")
    plt.scatter(yy, xx, color="green")
    plt.scatter(xx + yy, neg(yy + xx), color="green")
    plt.scatter(neg(xx + yy + xx + yy), yy + xx + neg(yy + xx), color="green")
    plt.plot([0, 1], [0, 1], "-k")
    plt.plot([0, 1], [0, 0], "-k")
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)

def plot_xz_grid(xyz):
    # project to x-z plane
    xy = list(set([(x, z) for x, y, z in xyz]))
    xx = [xi for (xi, yi) in xy]
    yy = [yi for (xi, yi) in xy]
    fig = plt.figure("porexz")
    fig.set_size_inches(4, 4)
    plt.scatter(neg(xx) + xx, yy + yy)

def plot_polygon(ax, polygon):
    settings = dict(closed=True, facecolor="#eeeeee", linewidth=1.,
                    edgecolor="black")
    polygon = np.array(polygon)
    polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])

    patch = patches.Polygon(polygon, **settings)
    patchm = patches.Polygon(polygon_m, **settings)
    #patch.set_zorder(10)
    #patchm.set_zorder(10)
    ax.add_patch(patch)
    ax.add_patch(patchm)

# will result in roughly nz * nr*(nr+1)/2 points
def tensorgrid(nz=30, nr=5, plot=False, eps=1e-2, eps2=8e-2, buf=10.,
               **params):
    params = Params(pugh_params) | Params(params)
    r = params.rMolecule
    r = r + eps

    # ---- create z part of tensor grid -----
    ztop = params.hpore/2.
    zbot = -ztop

    # 6 nodes => 5 sections
    z = [zbot - buf,
         zbot - r,
         ztop - params.h2 + r,
         ztop - params.h1 + r,
         ztop + r,
         ztop + buf]
    # relative meshwidths, radii
    hz = np.array([1., 1., .5, 1., 1.])
    rpore = [params.l0/2.,
             params.l3/2. - r,
             params.l2/2. - r,
             params.l1/2. - r,
             params.l0/2.]
    # to which of the two interval the shared endpoint belongs
    ep = [0, 1, 1, 1]

    grids = grid_piecewise1D(z, hz, N=nz, ep=ep)

    # ---- create xy (triangle) part of tensor grid -----
    # points in the unit triangle
    x = np.linspace(2*eps2, 1-eps, nr)
    y = np.linspace(eps2, 1-eps2, nr)
    xy = [(xi, yi) for xi in x for yi in y if xi > yi]

    # ---- tensor product
    xyz = tensor(xy, grids, rpore)

    if plot:
        print "Created %d points in z direction." % (sum(len(g) for g in grids),)
        print "Created %d points in xy direction." % (len(xy),)
        print "Total number of points:", len(xyz)
        #plot_1Dgrid(z, grids)
        plot_2Dgrid(xy)
        plot_xz_grid(xyz)
        ax = plt.gca()
        from nanopores.models.pughpore import polygon
        plot_polygon(ax, polygon())
    return xyz

if __name__ == "__main__":
    xyz = tensorgrid(nz=30, nr=5, eps2=3e-2, plot=True)
    plt.show()
#........................R.............................
#                                                     .
#                                                     .
#              .........l0..........                  .
#              .                   .                  .
#              ._ _______________ _...............    .
#              |D|               |D|     .   .   .    .
#              |D|......l1.......|D|    h1   .   .    .
#              |D|_ ____l2_____ _|D|......   h2  .    .
#              |DDD|_ _______ _|DDD|..........   .    .
#              |DDDDD|       |DDDDD|             .    .
#              |DDDDD|       |DDDDD|             .    .
#       DNA--->|DDDDD|       |DDDDD|           hpore  .
#              |DDDDD|       |DDDDD|             .    .
#              |DDDDD|..l3...|DDDDD|             .    .
#   MEMBRANE   |DDDDD|       |DDDDD|             .    H
#      |       |DDDDD|       |DDDDD|             .    .
#      |       |DDDDD|       |DDDDD|....h4       .    .
#______V_________|DDD|       |DDD|_____.________ .___ .......
#MMMMMMMMMMMMMMMM|DDD|       |DDD|MMMMM.MMMMMMMMM.MMMM.    hmem
#MMMMMMMMMMMMMMMM|DDD|_______|DDD|MMMMM.MMMMMMMMM.MMMM.......
#                .               .                    .
#                .......l4........                    .
#                                                     .
#                                                     .
#                                                     .
#......................................................