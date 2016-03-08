""" define solvers that record cumulative times needed for every loop """

from dolfin import *
from nanopores.tools.pdesystem import newtonsolve

def newton_solve(self, tol=None, damp=None, verbose=True):
    if tol is None: tol = self.tolnewton
    if damp is None: damp = self.newtondamp
    S = self.solvers.values()[0]    
    S.newtondamp = damp
    tcum = 0.
    Uold = self.solutions(deepcopy=True)
    
    for i in range(self.imax):
        tloop = Timer("loop")
        S.solve()
        if verbose:
            print '     Relative L2 Newton error:',S.relerror()
        if S.convergence(tol):
            if verbose:
                print "     Break loop because tolerance %s was reached." %tol
            converged = True
            break
        # cumulative time 
        tcum += tloop.stop()
        # calculate the error
        U = self.solutions(deepcopy=True)
        err = sum(errornorm(u, uold, "L2", degree_rise=0) for u, uold in zip(U, Uold)) / sum(norm(u, "L2") for u in U) 
        Uold = U
        self.save_estimate("err newton i", err, N=i+1)
        self.save_estimate("err newton time", err, N=tcum)
        continue
    else:
        if verbose: print "     Did not reach tolerance %s." %tol
        converged = False
    print "     Newton iterations:",i+1
        #print '     Relative L2 Newton error:',S.relerror()
    return i+1, converged
    

def hybrid_solve(self, tol=None, damp=None):
    if tol is None: tol = self.params["tolnewton"]
    if damp is None: damp = self.params["damp"]
    I = self.params["ipicard"]
    J = self.params["inewton"]
    nverbose = self.params["nverbose"]
    verbose = self.params["verbose"]
    times = {name : 0. for name in self.solvers}
    tcum = 0.
    Uold = self.solutions(deepcopy=True)

    for i in range(1, I+1):
        #v = self.functions["npm"]
        #plot(v)
        #self.visualize()
        if verbose:
            print "\n-- Fixed-Point Loop %d of max. %d" % (i, I)
        tloop = Timer("loop")
        for name, solver in self.solvers.items():
            if verbose:
                print "    Solving %s." % name
            t = Timer(name)
            if solver.is_linear:
                solver.solve()
                times[name] += t.stop()
            else:
                j, con = newtonsolve(solver, tol, damp, J, nverbose)
                times[name] += t.stop()
                if j==1 and con:
                    print "- Break at iteration %d because Newton stopped changing." %i
                    break
        else:
            # cumulative time 
            tcum += tloop.stop()
            # calculate the error
            U = self.solutions(deepcopy=True)
            err = sum(errornorm(u, uold, "L2", degree_rise=0) for u, uold in zip(U, Uold)) / sum(norm(u, "L2") for u in U) 
            Uold = U
            self.save_estimate("err hybrid i", err, N=i)
            self.save_estimate("err hybrid time", err, N=tcum)
            continue
        break
        
    Tt = sum(times.values())
    if verbose:
        print "\n CPU Time (solve): %.2f s" % Tt
        for name in self.solvers:
            print "  -) %s: %.2f s" %(name, times[name])
            
            
