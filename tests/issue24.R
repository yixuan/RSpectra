library(RSpectra)
library(Matrix)

## Test that eigs_sym() supports the symmetric Matrix types:
## dsyMatrix (dense), dsCMatrix (sparse, column), dsRMatrix (sparse, row)

set.seed(123)
n <- 100L
k <- 5L

## A symmetric dense matrix and its symmetric Matrix representations
A <- matrix(rnorm(n * n), n)
A <- A + t(A)
dsC <- as(as(A, "dgCMatrix"), "symmetricMatrix")  # dsCMatrix
smat <- list(
    dsyMatrix = as(A, "dsyMatrix"),
    dsCMatrix = dsC,
    dsRMatrix = as(dsC, "RsparseMatrix")
)
stopifnot(inherits(smat$dsyMatrix, "dsyMatrix"))
stopifnot(inherits(smat$dsCMatrix, "dsCMatrix"))
stopifnot(inherits(smat$dsRMatrix, "dsRMatrix"))

## Reference eigenvalues from the plain matrix
ref <- eigs_sym(A, k)

## Residual check: A %*% v ~= lambda * v
eresid <- function(eig)
    max(abs(A %*% eig$vectors - sweep(eig$vectors, 2, eig$values, "*")))

## eigs_sym() using the stored triangle (no conflict with 'uplo')
for (nm in names(smat)) {
    M <- smat[[nm]]
    res <- eigs_sym(M, k, lower = (M@uplo == "L"))
    stopifnot(max(abs(ref$values - res$values)) < 1e-8)
    stopifnot(eresid(res) < 1e-7)
}
cat("OK: eigs_sym() produces correct results for dsyMatrix, dsCMatrix, and dsRMatrix\n")

## A conflicting 'lower' should warn, and the 'uplo' slot must be used instead
warn <- function(expr) {
    w <- NULL
    withCallingHandlers(
        expr,
        warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") }
    )
    w
}
for (nm in names(smat)) {
    M <- smat[[nm]]
    ## Pass the triangle opposite to the one actually stored
    conf <- eigs_sym(M, k, lower = (M@uplo == "U"))
    stopifnot(!is.null(warn(eigs_sym(M, k, lower = (M@uplo == "U")))))
    ## Despite the warning, results still match the reference (uplo was used)
    stopifnot(max(abs(ref$values - conf$values)) < 1e-8)
}
cat("OK: eigs_sym() warns on 'lower'/'uplo' conflict and uses 'uplo'\n")
