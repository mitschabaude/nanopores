import dolfin
import nanopores
import colormaps as cm
import matplotlib
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

def porestreamlines(polygon=None, rx=10., ry=10., Nx=100, Ny=100, maxvalue=None, **fields):
    "streamlines plot of vector field around nanopore"  
    
    # interpolate on regular mesh symmetric w.r.t. center axis
    mesh2D = nanopores.RectangleMesh([-rx-0.1,-ry-0.1], [rx+0.1,ry+0.1], Nx, Ny)
    fields2 = nanopores.convert2D(mesh2D, *(fields.values()))
    
    # prepare polygon and copy to left half
    settings = dict(closed=True, facecolor="#eeeeee", linewidth=3., edgecolor="black")
    if polygon is not None:
        polygon = np.array(polygon)
        polygon_m = np.column_stack([-polygon[:,0], polygon[:,1]])

    # prepare plots
    Ny += 1
    Nx += 1
    Y, X = np.mgrid[-ry:ry:Ny*1j, -rx:rx:Nx*1j]
    U = np.zeros((Ny,Nx))
    V = np.zeros((Ny,Nx))
    formt = matplotlib.ticker.FuncFormatter(fmt)
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
        
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    if a=="1.0":
        return r"$10^{{{}}}$".format(b)
    elif a=="0.0":
        return r"$0$"
    else:
        return r'${}\cdot10^{{{}}}$'.format(a,b)
        
# do the plot with imported force field and polygon
from forcefield import F0, Fel0, Fdrag0, params, Howorka
print "parameters:", params

F, Fel, Fdrag = F0, Fel0, Fdrag0
poly = Howorka.polygon()

rx, ry = 6., 8.
porestreamlines(poly, rx, ry, Fel=Fel, Fdrag=Fdrag) #, maxvalue = 1e-11)

# modify plot for better output in paper
fig1 = plt.figure("Fel")
fig2 = plt.figure("Fdrag")
fig1.delaxes(fig1.axes[1])
fig2.axes[0].set_ylabel("")
fig2.axes[1].set_ylabel("force [N]")

# save to paper dir
from folders import NUMERICSFIGDIR as DIR
nanopores.savefigs("streamplot", DIR)

#plt.show()

