# (c) 2016 Gregor Mitscha-Baude
import dolfin
import nanopores
from . import colormaps as cm
import matplotlib
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def streamlines(polygon=None, patches=None, R=10., Htop=10., Hbot=10.,
                    Nx=100, Ny=100, figsize=(5, 5),
                    maxvalue=None, **fields):
    "streamlines plot of 2D vector field around nanopore"
    # TODO: make work for 3D vector field

    # interpolate on regular mesh symmetric w.r.t. center axis
    mesh2D = nanopores.RectangleMesh([-R-0.1,-Hbot-0.1], [R+0.1,Htop+0.1], Nx, Ny)
    fields2 = nanopores.convert2D(mesh2D, *(list(fields.values())))

    # prepare polygon and copy to left half
    settings = dict(closed=True, facecolor="#eeeeee", linewidth=3.,
                    edgecolor="black")
    if polygon is not None:
        polygon = np.array(polygon)
        polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])

    # prepare plots
    Ny += 1
    Nx += 1
    Y, X = np.mgrid[-Hbot:Htop:Ny*1j, -R:R:Nx*1j]
    U = np.zeros((Ny,Nx))
    V = np.zeros((Ny,Nx))
    formt = matplotlib.ticker.FuncFormatter(exp_format)
    ticks = [0] + [10**n for n in range(-15, -8)]

    # determine uniform color range from fields (maybe round to nearest 10-power)
    if maxvalue is None:
        maxvalue = max(dolfin.norm(F.vector(), "linf") for F in fields2)
        #maxvalue = 10**int(np.log10(maxvalue))

    for i, F in enumerate(fields2):
        Fstr = list(fields.keys())[i]
        fig, ax = plt.subplots(num=Fstr, figsize=figsize)

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
        strm = plt.streamplot(X, Y, U, V, arrowsize=1.0, linewidth=0.75, density=2.0,
                              cmap=cm.viridis, color=strength, norm=norm)
        #if i==len(fields2)-1:
        plt.colorbar(strm.lines, ticks=ticks, format=formt)
        fig.axes[1].set_ylabel("Force [N]")
        #plt.xlabel('x [nm]') #, fontsize=20)
        #plt.ylabel('z [nm]') #, fontsize=20)
        fig.axes[0].tick_params(
            axis="both",       # changes apply to both axes
            which="both",      # both major and minor ticks are affected
            bottom="off", top="off", left="off", right="off",
            labelbottom="off", labeltop="off",
            labelleft="off", labelright="off")
        

        # plot pore polygon on top
        if polygon is not None:
            patch = patches.Polygon(polygon, **settings)
            patchm = patches.Polygon(polygon_m, **settings)
            patch.set_zorder(10)
            patchm.set_zorder(10)
            ax.add_patch(patch)
            ax.add_patch(patchm)
            
        if patches is not None:
            for p in patches[i]:
                p.set_zorder(10)
                ax.add_patch(p)


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
    pass
