library(RSpectra)
library(Matrix)
n = 100
p = 50
k = 5

## Set up the test matrix
set.seed(123)
x = matrix(rnorm(n * p), n)
x[sample(n * p, floor(n * p / 2))] = 0

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

## Different matrix types
x1 = x
x2 = as(x, "dgeMatrix")
x3 = as(x, "dgCMatrix")
x4 = as(x, "dgRMatrix")

## Regular mode
svd0 = svd(x)
res1 = svds(x1, k)
res2 = svds(x2, k)
res3 = svds(x3, k)
res4 = svds(x4, k)
svd_resid(res1, svd0)
svd_resid(res2, svd0)
svd_resid(res3, svd0)
svd_resid(res4, svd0)

## Center x
xc = sweep(x, 2, colMeans(x), "-")
svd0 = svd(xc)
res1 = svds(x1, k, opts = list(center = TRUE))
res2 = svds(x2, k, opts = list(center = TRUE))
res3 = svds(x3, k, opts = list(center = TRUE))
res4 = svds(x4, k, opts = list(center = TRUE))
svd_resid(res1, svd0)
svd_resid(res2, svd0)
svd_resid(res3, svd0)
svd_resid(res4, svd0)

## Scale x
xs = sweep(x, 2, sqrt(colSums(x^2)), "/")
svd0 = svd(xs)
res1 = svds(x1, k, opts = list(scale = TRUE))
res2 = svds(x2, k, opts = list(scale = TRUE))
res3 = svds(x3, k, opts = list(scale = TRUE))
res4 = svds(x4, k, opts = list(scale = TRUE))
svd_resid(res1, svd0)
svd_resid(res2, svd0)
svd_resid(res3, svd0)
svd_resid(res4, svd0)

## Center and scale x
xcs = sweep(xc, 2, sqrt(colSums(xc^2)), "/")
svd0 = svd(xcs)
res1 = svds(x1, k, opts = list(center = TRUE, scale = TRUE))
res2 = svds(x2, k, opts = list(center = TRUE, scale = TRUE))
res3 = svds(x3, k, opts = list(center = TRUE, scale = TRUE))
res4 = svds(x4, k, opts = list(center = TRUE, scale = TRUE))
svd_resid(res1, svd0)
svd_resid(res2, svd0)
svd_resid(res3, svd0)
svd_resid(res4, svd0)

## Center and scale with given vectors
ctr = rnorm(p)
scl = abs(rnorm(p))
y = sweep(x, 2, ctr, "-")
y = sweep(y, 2, scl, "/")
svd0 = svd(y)
res1 = svds(x1, k, opts = list(center = ctr, scale = scl))
res2 = svds(x2, k, opts = list(center = ctr, scale = scl))
res3 = svds(x3, k, opts = list(center = ctr, scale = scl))
res4 = svds(x4, k, opts = list(center = ctr, scale = scl))
svd_resid(res1, svd0)
svd_resid(res2, svd0)
svd_resid(res3, svd0)
svd_resid(res4, svd0)
