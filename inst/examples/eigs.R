library(rARPACK)
library(Matrix)
n = 1000
k = 5

## Set up test matrices
set.seed(123)
x = matrix(rnorm(n^2), n)
x[sample(n^2, floor(n^2 / 2))] = 0
# General matrices
gen = list(x,
           as(x, "dgeMatrix"),
           as(x, "dgCMatrix"),
           as(x, "dgRMatrix"))
# Symmetric matrices
sym1 = list(x + t(x))
sym2 = list(as(x + t(x), "dsyMatrix"))

## Test whether the calculated eigenvalues and eigenvectors satisfy
##                         A * x = lambda * x
## Return the largest residual
eigen_resid = function(x, e)
{
    x = as.matrix(x)
    resid = x %*% e$vectors - e$vectors %*% diag(e$values)
    maxerr = max(abs(resid))
    return(paste("residual <", format(maxerr, digits = 5)))
}

## Capture test result, including error and warning
capture = function(expr, env)
{
    warn = NULL
    t1 = Sys.time()
    res = withCallingHandlers(
        tryCatch(eigen_resid(env$x, eval(expr, envir = env)),
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
gen_eigs_test = function(x, k)
{
    env = new.env()
    env$x = x
    env$k = k
    tests = expression(
        eigs(x, k, which = "LM"),
        eigs(x, k, which = "SM"),
        eigs(x, k, which = "LR"),
        eigs(x, k, which = "SR"),
        eigs(x, k, which = "LI"),
        eigs(x, k, which = "SI"),
        eigs(x, k, sigma = 0),
        eigs(x, k, sigma = 2),
        eigs(x, k, sigma = 1 + 1i)
    )
    res = lapply(tests, capture, env = env)
    cat(sprintf("[x of type '%s']:\n", class(x)))
    lapply(res, output)
    cat("\n")
    invisible(NULL)
}

invisible(lapply(gen, gen_eigs_test, k = k))

## Test symmetric matrices, eigs_sym() interface
sym1_eigs_test = function(x, k)
{
    env = new.env()
    env$x = x
    env$k = k
    tests = expression(
        eigs_sym(x, k, which = "LM"),
        eigs_sym(x, k, which = "SM"),
        eigs_sym(x, k, which = "LA"),
        eigs_sym(x, k, which = "SA"),
        eigs_sym(x, k, which = "BE"),
        eigs_sym(x, k, sigma = 0),
        eigs_sym(x, k, sigma = 2)
    )
    res = lapply(tests, capture, env = env)
    cat(sprintf("[x of type '%s']:\n", class(x)))
    lapply(res, output)
    cat("\n")
    invisible(NULL)
}

invisible(lapply(sym1, sym1_eigs_test, k = k))

## Test symmetric matrices, eigs() interface
sym2_eigs_test = function(x, k)
{
    env = new.env()
    env$x = x
    env$k = k
    tests = expression(
        eigs(x, k, which = "LM"),
        eigs(x, k, which = "SM"),
        eigs(x, k, which = "LA"),
        eigs(x, k, which = "SA"),
        eigs(x, k, which = "BE"),
        eigs(x, k, sigma = 0),
        eigs(x, k, sigma = 2)
    )
    res = lapply(tests, capture, env = env)
    cat(sprintf("[x of type '%s']:\n", class(x)))
    lapply(res, output)
    cat("\n")
    invisible(NULL)
}

invisible(lapply(sym2, sym2_eigs_test, k = k))
