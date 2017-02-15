from dolfin import has_lu_solver_method, has_krylov_solver_preconditioner

lusolver = "superlu_dist" if has_lu_solver_method("superlu_dist") else "default"

direct = dict(
    reuse = False,
    iterative = False,
    lusolver = lusolver,
)
direct_reuse = dict(
    reuse = True,
    iterative = False,
    lusolver = lusolver,
    luparams = dict(
        symmetric = False,
        same_nonzero_pattern = True,
        reuse_factorization = True,),
)

bicgstab = dict(
    reuse = False,
    iterative = True,
    lusolver = lusolver,
    luparams = dict(
        symmetric = False,
        same_nonzero_pattern = True,
        reuse_factorization = False,),
    ks = "bicgstab",
    kp = "hypre_euclid" if has_krylov_solver_preconditioner("hypre_euclid") else "ilu",
    kparams = dict(
        maximum_iterations = 600,
        monitor_convergence = False,
        nonzero_initial_guess = True,
        error_on_nonconvergence = True, #False,
        absolute_tolerance = 1e-5,
        relative_tolerance = 1e-8,
        preconditioner = dict(
            ilu = dict(fill_level = 1)))
)

poisson = dict(
    reuse = True,
    iterative = True,
    lusolver = lusolver,
    luparams = dict(
        symmetric = True,
        same_nonzero_pattern = True,
        reuse_factorization = True,),
    ks = "cg",
    kp = "hypre_euclid" if has_krylov_solver_preconditioner("hypre_euclid") else "ilu",
    kparams = dict(
        maximum_iterations = 500,
        monitor_convergence = False,
        relative_tolerance = 1e-6,
        preconditioner = dict(
            ilu = dict(fill_level = 1)))
)

stokes = dict(
    #reuse = False, # DEBUG
    reuse = True,
    iterative = False,
    lusolver = lusolver,
    luparams = dict(
        symmetric = True,
        same_nonzero_pattern = True,
        reuse_factorization = True,),
    ks = "tfqmr",
    kp = "hypre_euclid" if has_krylov_solver_preconditioner("hypre_euclid") else "ilu",
    fieldsplit = False, #True,
    kparams = dict(
        maximum_iterations = 1000,
        monitor_convergence = False,
        # large rel.tol. together with nonzero initial guess = bad idea!!!
        relative_tolerance = 1e-5,
        # absolute tolerance must not be too large compared with newton tol
        # (but also not too low since that would be inefficient)
        absolute_tolerance = 1e-5,
        nonzero_initial_guess = True,
        error_on_nonconvergence = False,
        preconditioner = dict(
            report = False,
            #structure = "same_nonzero_pattern",
            ilu = dict(fill_level = 1)))
)

gmres = dict(
    reuse = False,
    iterative = False,
    lusolver = lusolver,
    luparams = dict(
        symmetric = False,
        same_nonzero_pattern = True,
        reuse_factorization = False,),
    ks = "gmres",
    kp = "hypre_euclid",
    kparams = dict(
        maximum_iterations = 500,
        monitor_convergence = False,
        # large rel.tol. together with nonzero initial guess = bad idea!!!
        relative_tolerance = 1e-12,
        # absolute tolerance must not be too large compared with newton tol
        # (but also not too low since that would be inefficient)
        absolute_tolerance = 1e-5,
        nonzero_initial_guess = True,
        error_on_nonconvergence = False,
        preconditioner = dict(
            report = False,
            #structure = "same_nonzero_pattern",
            ilu = dict(fill_level = 1)))
)
