import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import dolfin

def fem2contour(u):
    mesh = u.function_space().mesh()
    v2d = dolfin.vertex_to_dof_map(u.function_space())
    
    # extract x and y coordinates of nodes
    x = mesh.coordinates()[:,0]
    y = mesh.coordinates()[:,1]
    triangles = mesh.cells()
    
    # Create triangulation.
    triang = mtri.Triangulation(x, y, triangles)
    
    # create array of node values from function
    z = u.vector()[v2d]

    # Plot the triangulation.
    plt.figure()
    plt.tricontourf(triang, z)
    #plt.triplot(triang, 'k-')
    #plt.title('Triangular grid')

if __name__ == "__main__":
    import HoworkaTools
    
    # get some mesh
    geo = HoworkaTools.geo
    mesh = geo.submesh("fluid")
    
    # Create some Function
    V = dolfin.FunctionSpace(mesh, "CG", 1)
    u = dolfin.Function(V)
    expr = dolfin.Expression("sin(x[0])")
    u.interpolate(expr)
    
    fem2contour(u)
    plt.show()
    
    
    
