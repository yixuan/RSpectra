## Internal function to test symmetry of matrices
## Intended to extend isSymmetric(), which has issues on dgRMatrix
is_sym <- function(mat)
{
    ## Early stop if nrow != ncol
    if (nrow(mat) != ncol(mat))
        return(FALSE)

    ## For dgCMatrix
    if (inherits(mat, "dgCMatrix"))
        return(.Call("is_sym_dgCMatrix", mat, 100 * .Machine$double.eps, PACKAGE = "RSpectra"))

    ## For dgRMatrix
    if (inherits(mat, "dgRMatrix"))
        return(.Call("is_sym_dgRMatrix", mat, 100 * .Machine$double.eps, PACKAGE = "RSpectra"))

    ## Default implementation
    isSymmetric(mat)
}



## Some simple tests
# library(Matrix)
# set.seed(123)
# x = rnorm(10000) * rbinom(10000, 1, 0.3)
# A = matrix(x, 100, 100)
# B = A + t(A)
#
# A1 = as.matrix(A)
# A2 = as(A, "dgCMatrix")
# A3 = as(A, "dgRMatrix")
# A4 = as(A, "dgeMatrix")
#
# B1 = as.matrix(B)
# B2 = as(B, "dgCMatrix")
# B3 = as(B, "dgRMatrix")
# B4 = as(B, "dgeMatrix")
# B5 = as(B, "dsyMatrix")
#
# RSpectra:::is_sym(A1)
# RSpectra:::is_sym(A2)
# RSpectra:::is_sym(A3)
# RSpectra:::is_sym(A4)
#
# RSpectra:::is_sym(B1)
# RSpectra:::is_sym(B2)
# RSpectra:::is_sym(B3)
# RSpectra:::is_sym(B4)
# RSpectra:::is_sym(B5)
