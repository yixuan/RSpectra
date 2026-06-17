library(RSpectra)
library(Matrix)

## Test that eigs_sym() supports the symmetric sparse types dsCMatrix and dsRMatrix

set.seed(123)
n <- 100L
k <- 5L

## A symmetric dense matrix and its symmetric sparse representations
A <- matrix(rnorm(n * n), n)
A <- A + t(A)
dsC <- as(as(A, "dgCMatrix"), "symmetricMatrix")  # dsCMatrix
dsR <- as(dsC, "RsparseMatrix")                   # dsRMatrix
stopifnot(inherits(dsC, "dsCMatrix"), inherits(dsR, "dsRMatrix"))

## Reference eigenvalues from the dense matrix
ref <- eigs_sym(A, k)

## eigs_sym() using the stored triangle (no conflict with 'uplo')
resC <- eigs_sym(dsC, k, lower = (dsC@uplo == "L"))
resR <- eigs_sym(dsR, k, lower = (dsR@uplo == "L"))

stopifnot(max(abs(ref$values - resC$values)) < 1e-8)
stopifnot(max(abs(ref$values - resR$values)) < 1e-8)

## Residual check: A %*% v ~= lambda * v
eresid <- function(eig)
    max(abs(A %*% eig$vectors - sweep(eig$vectors, 2, eig$values, "*")))
stopifnot(eresid(resC) < 1e-7)
stopifnot(eresid(resR) < 1e-7)
cat("OK: eigs_sym() produces correct results for dsCMatrix and dsRMatrix\n")

## A conflicting 'lower' should warn, and the 'uplo' slot must be used instead
warn <- function(expr) {
    w <- NULL
    withCallingHandlers(
        expr,
        warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") }
    )
    w
}

## Pass the triangle opposite to the one actually stored
confC <- eigs_sym(dsC, k, lower = (dsC@uplo == "U"))
confR <- eigs_sym(dsR, k, lower = (dsR@uplo == "U"))
stopifnot(!is.null(warn(eigs_sym(dsC, k, lower = (dsC@uplo == "U")))))
stopifnot(!is.null(warn(eigs_sym(dsR, k, lower = (dsR@uplo == "U")))))

## Despite the warning, results still match the reference (uplo was used)
stopifnot(max(abs(ref$values - confC$values)) < 1e-8)
stopifnot(max(abs(ref$values - confR$values)) < 1e-8)
cat("OK: eigs_sym() warns on 'lower'/'uplo' conflict and uses 'uplo'\n")
