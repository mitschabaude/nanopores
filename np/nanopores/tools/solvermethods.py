from dolfin import *

direct = dict(
    reuse = False,
    iterative = False,
    lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
)
direct_reuse = dict(
    reuse = True,
    iterative = False,
    lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
)


bicgstab = dict(
    reuse = False,
    iterative = True,
    lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
    ks = "bicgstab",
    kp = "hypre_euclid",
    kparams = dict(
        maximum_iterations = 500,
        monitor_convergence = False,
        relative_tolerance = 1e-6,
        preconditioner = dict(
            ilu = dict(fill_level = 1)))
)

poisson = dict(
    reuse = True,
    iterative = True,
    lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
    luparams = dict(
        symmetric = True,
        same_nonzero_pattern = True,
        reuse_factorization = True,),
    ks = "cg",
    kp = "hypre_euclid",
    kparams = dict(
        maximum_iterations = 500,
        monitor_convergence = False,
        relative_tolerance = 1e-6,
        preconditioner = dict(
            ilu = dict(fill_level = 1)))
)

stokes = dict(
    reuse = True,
    iterative = False,
    lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
    luparams = dict(
        symmetric = True,
        same_nonzero_pattern = True,
        reuse_factorization = True,),
    ks = "tfqmr",
    kp = "hypre_euclid",
    fieldsplit = False, #True,
    kparams = dict(
        maximum_iterations = 5000,
        monitor_convergence = True,
        # large rel.tol. together with nonzero initial guess = bad idea!!!
        relative_tolerance = 1e-8,
        # absolute tolerance must not be too large compared with newton tol
        # (but also not too low since that would be inefficient)
        absolute_tolerance = 1e-5,
        nonzero_initial_guess = True,
        error_on_nonconvergence = False,
        preconditioner = dict(
            report = False,
            structure = "same_nonzero_pattern",
            ilu = dict(fill_level = 1)))
)
