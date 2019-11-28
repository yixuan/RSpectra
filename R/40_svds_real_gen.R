svds_real_gen <- function(A, k, nu, nv, opts, mattype, extra_args = list())
{
    if (mattype == "function")
    {
        m = as.integer(extra_args$dim[1])
        n = as.integer(extra_args$dim[2])
    } else {
        m = nrow(A)
        n = ncol(A)
    }
    wd = min(m, n)

    # Check for matrices that are too small
    if (wd < 3)
        stop("nrow(A) and ncol(A) should be at least 3")

    # If all singular values are requested, call svd() instead,
    # and give a warning
    if (k == wd)
    {
        warning("all singular values are requested, svd() is used instead")
        return(c(svd(A, nu = nu, nv = nv),
               nconv = wd, niter = 0))
    }

    # Matrix will be passed to C++, so we need to check the type.
    # Convert the matrix type if A is stored other than double.
    #
    # However, for sparse matrices defined in Matrix package,
    # they are always double, so we can omit this check.
    if (mattype == "matrix" & typeof(A) != "double")
    {
        mode(A) = "double"
    }

    # Check the value of 'k'
    if (k <= 0 | k >= wd)
        stop("'k' must satisfy 0 < k < min(nrow(A), ncol(A)).\nTo calculate all singular values, try svd()")

    # Check the values of 'nu' and 'nv'
    if (nu < 0 | nv < 0 | nu > k | nv > k)
        stop("'nu' and 'nv' must satisfy 0 <= nu <= k and 0 <= nv <= k")

    # Arguments to be passed to Spectra
    spectra.param = list(ncv = min(wd, max(2 * k + 1, 20)),
                         tol = 1e-10,
                         maxitr = 1000,
                         center = FALSE,
                         scale = FALSE)
    # By default center = FALSE and scale = FALSE
    ctr = rep(0, n)
    scl = rep(1, n)

    # Update ctr and scl from opts
    # 1. If `center == TRUE`, then the centering vector is the column mean of A
    # 2. If `center` is a vector, then use this vector to center A
    # 3. In other cases, do not center A
    if (isTRUE(opts$center))
    {
        ctr = colMeans(A)
    } else if (is.numeric(opts$center)) {
        if (length(opts$center) != n)
            stop("opts$center must be TRUE/FALSE or a vector of length n")

        ctr = as.numeric(opts$center)
        opts$center = TRUE
    } else {
        opts$center = FALSE
    }
    # Scaling is always applied to vectors **after centering**
    # 4. If `scale == TRUE`, then the scaling vector consists of the norms of column
    #    vectors of A **after centering**
    # 5. If `scale` is a vector, then use this vector to scale A
    # 6. In other cases, do not scale A
    if (isTRUE(opts$scale))
    {
        sumx = colSums(A)
        sumxx = colSums(A^2)
        scl = sqrt(sumxx - 2 * sumx * ctr + m * ctr^2)
    } else if (is.numeric(opts$scale)) {
        if (length(opts$scale) != n)
            stop("opts$scale must be TRUE/FALSE or a vector of length n")

        scl = as.numeric(opts$scale)
        opts$scale = TRUE
    } else {
        opts$scale = FALSE
    }

    # Update parameters from 'opts' argument
    spectra.param[names(opts)] = opts

    # Any other arguments passed to C++ code
    spectra.param = c(spectra.param, as.list(extra_args),
                      list(ctr_vec = ctr, scl_vec = scl))

    # Check the value of 'ncv'
    if (spectra.param$ncv <= k | spectra.param$ncv > wd)
        stop("'opts$ncv' must be > k and <= min(nrow(A), ncol(A))")

    # Call the C++ function
    res = .Call("svds_gen",
                A,
                as.integer(m), as.integer(n),
                as.integer(k), as.integer(nu), as.integer(nv),
                as.list(spectra.param),
                as.integer(MAT_TYPE[mattype]),
                PACKAGE = "RSpectra")

    return(res)
}
