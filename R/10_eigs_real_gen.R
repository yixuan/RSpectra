eigs.real_gen <- function(A, n, k, which, sigma, opts, ..., mattype,
                          extra_args = list())
{
    # Check whether 'A' is a square matrix
    # Skip this step if A is a function
    if (!is.null(dim(A)))
    {
        if (nrow(A) != ncol(A) | nrow(A) != n)
            stop("'A' must be a square matrix of size n")
    }
    
    # eigs() is not suitable for small matrices
    if (n < 3)
        stop("dimension of 'A' must be at least 3")

    # If >= n-1 eigenvalues are requested, call eigen() instead,
    # and give a warning
    if (k == n - 1)
    {
        res = eigen(A, only.values = identical(opts$retvec, FALSE))
        exclude = switch(which,
                         LM = n,
                         SM = 1,
                         LR = which.min(Re(res$values)),
                         SR = which.max(Re(res$values)),
                         LI = which.min(abs(Im(res$values))),
                         SI = which.max(abs(Im(res$values))))
        warning("(n-1) eigenvalues are requested, eigen() is used instead")
        return(list(values = res$values[-exclude],
                    vectors = res$vectors[, -exclude],
                    nconv = n, niter = 0))
    }
    if (k == n)
    {
        res = eigen(A, only.values = identical(opts$retvec, FALSE))
        warning("all eigenvalues are requested, eigen() is used instead")
        return(c(res, nconv = n, niter = 0))
    }

    # Matrix will be passed to C++, so we need to check the type.
    # ARPACK only supports matrices in float or double, so we need
    # to do the conversion if A is stored other than double.
    #
    # However, for sparse matrices defined in Matrix package,
    # they are always double, so we can omit this check.
    if (mattype == "matrix" & typeof(A) != "double")
    {
        mode(A) = "double"
    }
    # Check the value of 'k'
    if (k <= 0 | k >= n - 1)
        stop("'k' must satisfy 0 < k < nrow(A) - 1.\nTo calculate all eigenvalues, try eigen()")

    # Check sigma
    # workmode == 1: ordinary
    # workmode == 3: Shift-invert mode
    if (is.null(sigma))
    {
        workmode = "regular"
        sigma = 0
    } else {
        if(abs(Im(sigma)) < 1e-16)
        {
            sigma = Re(sigma)
            workmode = "real_shift"
        } else {
            workmode = "complex_shift"
        }
    }

    # Arguments to be passed to ARPACK
    arpack.param = list(which = which,
                        ncv = min(n, max(2 * k + 1, 20)),
                        tol = 1e-10,
                        maxitr = 1000,
                        retvec = TRUE,
                        sigmar = Re(sigma[1]),
                        sigmai = Im(sigma[1]))

    # Check the value of 'which'
    eigenv.type = c("LM", "SM", "LR", "SR", "LI", "SI")
    if (!(arpack.param$which %in% eigenv.type))
    {
        stop(sprintf("argument 'which' must be one of\n%s",
                     paste(eigenv.type, collapse = ", ")))
    }

    # Update parameters from 'opts' argument
    arpack.param[names(opts)] = opts
    arpack.param$which = EIGS_RULE[arpack.param$which]
    
    # Any other arguments passed to C++ code, for example use_lower and fun_args
    arpack.param = c(arpack.param, as.list(extra_args))

    # Check the value of 'ncv'
    if (arpack.param$ncv < k + 2 | arpack.param$ncv > n)
        stop("'opts$ncv' must be >= k+2 and <= nrow(A)")

    # Call the C++ function
    fun = switch(workmode,
                 regular = "eigs_gen",
                 real_shift = "eigs_real_shift_gen",
                 complex_shift = "eigs_complex_shift_gen",
                 stop("unknown work mode"))
    res = .Call(fun,
                A,
                as.integer(n), as.integer(k),
                as.list(arpack.param),
                as.integer(MAT_TYPE[mattype]),
                PACKAGE = "rARPACK")

    return(res)
}
