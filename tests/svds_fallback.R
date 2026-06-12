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

## Test that svds() with function interface works when k == min(m, n)
## (previously the svd() fallback would error because A is a closure)

set.seed(123)
m <- 8
n <- 5
k <- min(m, n)
A <- matrix(rnorm(m * n), m)
Atrans <- t(A)
s0 <- svd(A)

fun_A  <- function(x, args) as.numeric(args %*% x)
fun_At <- function(x, args) as.numeric(crossprod(args, x))

s1 <- svds(fun_A, k, Atrans = fun_At, dim = c(m, n), args = A)

stopifnot(max(abs(s0$d - s1$d)) < 1e-10)
cat("OK: svds() function interface matches svd() when k == min(m, n)\n")

## Test that svds() with function interface and center/scale works when k == min(m, n)

set.seed(456)
m <- 8
n <- 5
k <- min(m, n)
A <- matrix(rnorm(m * n), m)
ctr <- runif(n)
scl <- runif(n, min = 0.5, max = 2)

fun_A  <- function(x, args) as.numeric(args %*% x)
fun_At <- function(x, args) as.numeric(crossprod(args, x))

s0 <- svd(sweep(sweep(A, 2, ctr), 2, scl, "/"))
s1 <- svds(fun_A, k, Atrans = fun_At, dim = c(m, n), args = A,
           opts = list(center = ctr, scale = scl))

stopifnot(max(abs(s0$d - s1$d)) < 1e-10)
cat("OK: svds() function interface with center/scale matches svd() when k == min(m, n)\n")
