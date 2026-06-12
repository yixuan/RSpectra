library(RSpectra)
library(Matrix)

## Test that svds() with center/scale gives correct result when k == min(m, n)
## (previously the svd() fallback ignored center and scale)

set.seed(42)
m <- 10
n <- 5
k <- min(m, n)
A <- matrix(rnorm(m * n), m)
ctr <- runif(n)
scl <- runif(n, min = 0.5, max = 2)

s0 <- svd(sweep(sweep(A, 2, ctr), 2, scl, "/"))
s1 <- svds(A, k, opts = list(center = ctr, scale = scl))

stopifnot(max(abs(s0$d - s1$d)) < 1e-10)
cat("OK: svds() with center/scale matches svd() when k == min(m, n)\n")
