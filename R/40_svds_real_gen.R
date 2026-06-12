# Fallback for full SVD when A is a function interface (k == min(m, n)).
# Recover A'A or AA' by applying A/Atrans to standard basis vectors,
# then use eigen(symmetric=TRUE) to obtain singular values and vectors.
svds_fallback_A_Atrans <- function(A, Atrans, m, n, nu, nv, fun_args)
{
    wd = min(m, n)
    nu = min(nu, wd)
    nv = min(nv, wd)

    if (m > n)
    {
        # Form A'A (n x n): column j = Atrans(A(e_j))
        AtA = matrix(0, n, n)
        for (j in seq_len(n))
        {
            ej = numeric(n)
            ej[j] = 1
            AtA[, j] = Atrans(A(ej, fun_args), fun_args)
        }
        # Eigenvalues are sigma^2, eigenvectors are V
        eig = eigen(AtA, symmetric = TRUE)
        d = sqrt(pmax(eig$values, 0))
        V = eig$vectors
        # Compute U = A * V * diag(1/d)
        # Only the first nu columns are needed
        if (nu > 0)
        {
            U = matrix(0, m, nu)
            for (i in seq_len(nu))
            {
                if (d[i] > 0)
                    U[, i] = A(V[, i], fun_args) / d[i]
            }
        } else {
            U = NULL
        }
    } else {
        # Form AA' (m x m): column j = A(Atrans(e_j))
        AAt = matrix(0, m, m)
        for (j in seq_len(m))
        {
            ej = numeric(m)
            ej[j] = 1
            AAt[, j] = A(Atrans(ej, fun_args), fun_args)
        }
        eig = eigen(AAt, symmetric = TRUE)
        d = sqrt(pmax(eig$values, 0))
        U = eig$vectors
        # Compute V = Atrans * U * diag(1/d)
        # Only the first nv columns are needed
        if (nv > 0)
        {
            V = matrix(0, n, nv)
            for (i in seq_len(nv))
            {
                if (d[i] > 0)
                    V[, i] = Atrans(U[, i], fun_args) / d[i]
            }
        } else {
            V = NULL
        }
    }

    list(d = d,
         u = if (nu > 0) U else NULL,
         v = if (nv > 0) V else NULL,
         nconv = wd,
         niter = 0)
}

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

    # If all singular values are requested, call svd() instead,
    # and give a warning
    if (k == wd)
    {
        warning("all singular values are requested, svd() is used instead")
        if (mattype == "function") {
            return(svds_fallback_A_Atrans(A, extra_args$Atrans,
                                          m, n, nu, nv,
                                          extra_args$fun_args))
        }
        # Apply centering and scaling if requested: B = (A - 1c')S
        Asvds = A
        if (isTRUE(opts$center))
            Asvds = sweep(Asvds, 2, ctr)
        if (isTRUE(opts$scale))
            Asvds = sweep(Asvds, 2, scl, "/")
        return(c(svd(Asvds, nu = nu, nv = nv),
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
