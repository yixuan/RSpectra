eigs.real_sym <- function(A, n, k, which, sigma, opts, ..., mattype,
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
    
    # If all eigenvalues are requested, call eigen() instead,
    # and give a warning
    if (k == n)
    {
        warning("all eigenvalues are requested, eigen() is used instead")
        return(c(eigen(if(extra_args$use_lower) A else t(A),
                       symmetric = TRUE,
                       only.values = identical(opts$retvec, FALSE)),
                 nconv = n, niter = 0))
    }
    
    # Matrix will be passed to C++, so we need to check the type.
    # ARPACK only supports matrices in float or double, so we need
    # to do the conversion if A is stored other than double.
    #
    # However, for dsyMatrix matrices defined in Matrix package,
    # they are always double, so we can omit this check. 
    if (mattype == "matrix" & typeof(A) != "double")
    {
        mode(A) = "double"
    }
    # Check the value of 'k'
    if (k <= 0 | k >= n)
        stop("'k' must satisfy 0 < k < nrow(A)")
    
    # Check sigma
    if (is.null(sigma))
    {
        workmode = "regular"
        sigma = 0
    } else {
        workmode = "real_shift"
        if(is.complex(sigma)) warning("only real part of sigma is used")
        sigma = Re(sigma)
    }
    
    # Arguments to be passed to ARPACK
    arpack.param = list(which = which,
                        ncv = min(n, max(2 * k + 1, 20)),
                        tol = 1e-10,
                        maxitr = 1000,
                        retvec = TRUE,
                        sigma = sigma)
    
    # Check the value of 'which'
    eigenv.type = c("LM", "SM", "LA", "SA", "BE")
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
    if (arpack.param$ncv <= k | arpack.param$ncv > n)
        stop("'opts$ncv' must be > k and <= nrow(A)")
    
    # Call the C++ function
    fun = switch(workmode,
                 regular = "eigs_sym",
                 real_shift = "eigs_shift_sym",
                 stop("unknown work mode"))
    res = .Call(fun,
                A,
                as.integer(n), as.integer(k),
                as.list(arpack.param),
                as.integer(MAT_TYPE[mattype]),
                PACKAGE = "rARPACK")
    
    return(res)
}
