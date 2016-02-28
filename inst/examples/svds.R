library(rARPACK)
library(Matrix)
n = 1000
p = 500
k = 5

## Set up test matrices
set.seed(123)
x = matrix(rnorm(n * p), n)
x[sample(n * p, floor(n * p / 2))] = 0
# General matrices
gen = list(x,
           as(x, "dgeMatrix"),
           as(x, "dgCMatrix"),
           as(x, "dgRMatrix"))
gent = lapply(gen, t)
# Symmetric matrices
sym = list(as(crossprod(x), "dsyMatrix"))

## Test whether the calculated (d, u, v) are consistent with svd()
## Return the largest residual
svd_resid = function(res, svd0)
{
    d_resid = svd0$d[1:length(res$d)] - res$d
    u_resid = v_resid = 0
    if(!is.null(res$u))
        u_resid = abs(svd0$u[, 1:ncol(res$u)]) - abs(res$u)
    if(!is.null(res$v))
        v_resid = abs(svd0$v[, 1:ncol(res$v)]) - abs(res$v)
    mabs = function(x) max(abs(x))
    maxerr = max(mabs(d_resid), mabs(u_resid), mabs(v_resid))
    return(paste("residual <", format(maxerr, digits = 5)))
}
# "True" values
gen0 = svd(x)
gen0t = svd(t(x))
sym0 = svd(crossprod(x))

## Capture test result, including error and warning
capture = function(expr, env)
{
    warn = NULL
    t1 = Sys.time()
    res = withCallingHandlers(
        tryCatch(svd_resid(eval(expr, envir = env), env$svd0),
                 error = function(e) e),
        warning = function(w) {warn <<- w; invokeRestart("muffleWarning")}
    )
    t2 = Sys.time()
    return(list(src = deparse(expr), res = res, warn = warn,
                time = as.numeric(t2 - t1)))
}

## Output result
output = function(res)
{
    cat(res$src, rep(" ", 32 - nchar(res$src)), ": ", sep = "")
    if(inherits(res$res, "error"))
    {
        cat("ERROR")
    } else cat(res$res)
    
    if(!is.null(res$warn)) cat(" (with warning)")
    cat("\n", rep(" ", 34), sprintf("(%f seconds)\n", res$time), sep = "")
}



## Test general matrices
svds_test = function(x, k, svd0)
{
    env = new.env()
    env$x = x
    env$k = k
    env$svd0 = svd0
    tests = expression(
        svds(x, k),
        svds(x, k, nu = 0),
        svds(x, k, nv = 0),
        svds(x, k, nu = 0, nv = 0)
    )
    res = lapply(tests, capture, env = env)
    cat(sprintf("[x of type '%s']:\n", class(x)))
    lapply(res, output)
    cat("\n")
    invisible(NULL)
}

invisible(lapply(gen, svds_test, k = k, svd0 = gen0))
invisible(lapply(gent, svds_test, k = k, svd0 = gen0t))
invisible(lapply(sym, svds_test, k = k, svd0 = sym0))
