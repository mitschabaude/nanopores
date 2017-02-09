# (c) 2016 Gregor Mitscha-Baude
import dolfin
import nanopores
import colormaps as cm
import matplotlib
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def streamlines(polygon=None, rx=10., ry=10., Nx=100, Ny=100,
                    maxvalue=None, **fields):
    "streamlines plot of 2D vector field around nanopore"
    # TODO: make work for 3D vector field

    # interpolate on regular mesh symmetric w.r.t. center axis
    mesh2D = nanopores.RectangleMesh([-rx-0.1,-ry-0.1], [rx+0.1,ry+0.1], Nx, Ny)
    fields2 = nanopores.convert2D(mesh2D, *(fields.values()))

    # prepare polygon and copy to left half
    settings = dict(closed=True, facecolor="#eeeeee", linewidth=3.,
                    edgecolor="black")
    if polygon is not None:
        polygon = np.array(polygon)
        polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])

    # prepare plots
    Ny += 1
    Nx += 1
    Y, X = np.mgrid[-ry:ry:Ny*1j, -rx:rx:Nx*1j]
    U = np.zeros((Ny,Nx))
    V = np.zeros((Ny,Nx))
    formt = matplotlib.ticker.FuncFormatter(exp_format)
    ticks = [0] + [10**n for n in range(-16, -9)]

    # determine uniform color range from fields (maybe round to nearest 10-power)
    if maxvalue is None:
        maxvalue = max(dolfin.norm(F.vector(), "linf") for F in fields2)
        #maxvalue = 10**int(np.log10(maxvalue))

    for i, F in enumerate(fields2):
        Fstr = fields.keys()[i]
        fig, ax = plt.subplots(figsize=(5, 4.5), num=Fstr)

        # fill array with function values
        for y in range(Ny):
            for x in range(Nx):
                f = F(X[y][x], Y[y][x])
                U[y][x] = f[0]
                V[y][x] = f[1]

        # streamplot with logarithmic scale
        strength = np.sqrt(U*U+V*V)
        norm = matplotlib.colors.SymLogNorm(linthresh=ticks[1], linscale=1.0,
                                            vmin=0., vmax=maxvalue)
        strm = plt.streamplot(X, Y, U, V, arrowsize=1.5, linewidth=1.5, density=1.5,
                              cmap=cm.viridis, color=strength, norm=norm)
        plt.colorbar(strm.lines, ticks=ticks, format=formt)
        plt.xlabel('x [nm]') #, fontsize=20)
        plt.ylabel('z [nm]') #, fontsize=20)

        # plot pore polygon on top
        if polygon is not None:
            patch = patches.Polygon(polygon, **settings)
            patchm = patches.Polygon(polygon_m, **settings)
            patch.set_zorder(10)
            patchm.set_zorder(10)
            ax.add_patch(patch)
            ax.add_patch(patchm)

def exp_format(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if a=="1.0":
        return r"$10^{{{}}}$".format(b)
    elif a=="0.0":
        return r"$0$"
    else:
        return r'${}\cdot10^{{{}}}$'.format(a,b)

if __name__ == "__main__":
    import nanopores.models.pughpore as pugh
    from nanopores.models.pughpoints import plot_polygon
    from nanopores.models.diffusion import get_pugh_diffusivity

    from nanopores.tools.utilities import uCross, RectangleMesh
    from math import pi, sqrt

    dparams = {2: dict(diamPore=6.*sqrt(pi)/2., diamDNA=2.5*sqrt(pi)/2.,
                       Nmax=1.2e5, dim=2, r=0.11, h=.75,
                       cheapest=False, Membraneqs=-.5),
               3: dict(diamPore=6., Nmax=1e6, dim=3, r=0.11, h=2.0, cheapest=False)}

    # obtain diffusivity field and project to x-z plane
    functions = get_pugh_diffusivity(**dparams[2])
    D3D = functions["D"][0]
    D0 = nanopores.D
    def F(x, z):
        if x>=0:
            return D3D([x, z])/D0
        else:
            return D3D([-x, z])/D0
    #D = uCross(u=D3D, axis=1, degree=1, dim=2)

    # obtain 2D mesh where we will evaluate field
    rx, ry = pugh.pughpore.params["R"], 0.5*pugh.pughpore.params["H"]
    rx, ry = 13, 28
    Nx, Ny = 101, 201
    #mesh2D = RectangleMesh([-R,-H/2.], [R, H/2.], int(4*R), int(2*H))

    Y, X = np.mgrid[-ry:ry:Ny*1j, -rx:rx:Nx*1j]
    U = np.zeros((Ny,Nx))

    for y in range(Ny):
        for x in range(Nx):
            U[y][x] = F(X[y][x], Y[y][x])

    fig, ax = plt.subplots(figsize=(6.5, 6.5), num="r0.11")
    pc = plt.pcolor(X, Y, U, cmap=cm.inferno, vmin=0, vmax=1)
    plt.colorbar(pc)
    plot_polygon(ax, pugh.polygon(diamPore=6., rmem=13))
    plt.xlim(-13, 13)
    plt.ylim(-25, 28)
    #plt.xlabel("x [nm]")
    #plt.ylabel("z [nm]")
    fig.axes[1].set_ylabel(r"$D_x / D_0$")
    #cb = fig.colorbar(CS, cax=cax, extend="both", orientation="vertical", format=formt)
    import folders
    nanopores.savefigs("pugh_Dfield", folders.FIGDIR)
    #plt.show()

