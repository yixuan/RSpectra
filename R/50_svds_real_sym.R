svds.real_sym <- function(A, k, nu, nv, opts, ..., mattype,
                          extra_args = list())
{
    n = nrow(A)
 
    # Check for matrices that are too small
    if (n < 3)
        stop("nrow(A) and ncol(A) should be at least 3")
    
    # If all singular values are requested, call svd() instead,
    # and give a warning
    if (k == n)
    {
        warning("all singular values are requested, svd() is used instead")
        return(c(svd(A, nu = nu, nv = nv),
                 nconv = n, niter = 0))
    }
    
    # Matrix will be passed to C++, so we need to check the type.
    # ARPACK only supports matrices in float or double, so we need
    # to do the conversion if A is stored other than double.
    #
    # However, for symmetric matrices defined in Matrix package,
    # they are always double, so we can omit this check. 
    if (mattype == "matrix" & typeof(A) != "double")
    {
        mode(A) = "double"
    }
    
    # Check the value of 'k'
    if (k <= 0 | k >= n)
        stop("'k' must satisfy 0 < k < min(nrow(A), ncol(A)).\nTo calculate all singular values, try svd()")
    
    # Check the values of 'nu' and 'nv'
    if (nu < 0 | nv < 0 | nu > k | nv > k)
        stop("'nu' and 'nv' must satisfy 0 <= nu <= k and 0 <= nv <= k")
    
    # Arguments to be passed to ARPACK
    arpack.param = list(ncv = min(n, max(2 * k + 1, 20)),
                        tol = 1e-10,
                        maxitr = 1000)
    
    # Update parameters from 'opts' argument
    arpack.param[names(opts)] = opts
    
    # Any other arguments passed to C++ code, for example use_lower
    arpack.param = c(arpack.param, as.list(extra_args))
    
    # Check the value of 'ncv'
    if (arpack.param$ncv <= k | arpack.param$ncv > n)
        stop("'opts$ncv' must be > k and <= min(nrow(A), ncol(A))")
    
    # Call the C++ function
    res = .Call("svds_sym",
                A,
                as.integer(n),
                as.integer(k), as.integer(nu), as.integer(nv),
                as.list(arpack.param),
                as.integer(MAT_TYPE[mattype]),
                PACKAGE = "rARPACK")
    
    return(res)
}
