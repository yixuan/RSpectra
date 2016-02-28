##' Find the Largest k Singular Values/Vectors of a Matrix
##'
##' @description
##' Given an \eqn{m} by \eqn{n} matrix \eqn{A},
##' function \code{svds()} can find its largest \eqn{k}
##' singular values and the corresponding singular vectors.
##' It is also called the Truncated Singular Value Decomposition
##' since it only contains a subset of the whole singular triplets.
##' 
##' Currently \code{svds()} supports matrices of the following classes:
##' 
##' \tabular{ll}{
##'   \code{matrix}     \tab The most commonly used matrix type,
##'                          defined in \strong{base} package.\cr
##'   \code{dgeMatrix}  \tab General matrix, equivalent to \code{matrix},
##'                          defined in \strong{Matrix} package.\cr
##'   \code{dgCMatrix}  \tab Column oriented sparse matrix, defined in
##'                          \strong{Matrix} package.\cr
##'   \code{dgRMatrix}  \tab Row oriented sparse matrix, defined in
##'                          \strong{Matrix} package.\cr
##'   \code{dsyMatrix}  \tab Symmetrix matrix, defined in \strong{Matrix}
##'                          package.
##' }
##'
##' Note that when \eqn{A} is symmetric,
##' SVD reduces to eigen decomposition, so you may consider using
##' \code{\link{eigs}()} instead.
##' 
##' @param A The matrix whose truncated SVD is to be computed.
##' @param k Number of singular values requested.
##' @param nu Number of left singular vectors to be computed. This must
##'           be between 0 and \code{k}.
##' @param nv Number of right singular vectors to be computed. This must
##'           be between 0 and \code{k}.
##' @param opts Control parameters related to the computing
##'             algorithm. See \strong{Details} below.
##' @param \dots Currently not used.
##'
##' @details The \code{opts} argument is a list that can supply any of the
##' following parameters:
##'
##' \describe{
##' \item{\code{ncv}}{Number of Lanzcos basis vectors to use. More vectors
##'                   will result in faster convergence, but with greater
##'                   memory use. \code{ncv} must be satisfy
##'                   \eqn{k < ncv \le p}{k < ncv <= p} where
##'                   \code{p = min(m, n)}.
##'                   Default is \code{min(p, max(2*k+1, 20))}.}
##' \item{\code{tol}}{Precision parameter. Default is 1e-10.}
##' \item{\code{maxitr}}{Maximum number of iterations. Default is 1000.}
##' }
##' 
##' @return A list with the following components:
##' \item{d}{A vector of the computed singular values.}
##' \item{u}{An \code{m} by \code{nu} matrix whose columns contain
##'          the left singular vectors. If \code{nu == 0}, \code{NULL}
##'          will be returned.}
##' \item{v}{An \code{n} by \code{nv} matrix whose columns contain
##'          the right singular vectors. If \code{nv == 0}, \code{NULL}
##'          will be returned.}
##' \item{nconv}{Number of converged singular values.}
##' \item{niter}{Number of iterations used.}
##' \item{nops}{Number of matrix-vector multiplications used.}
##' @author Yixuan Qiu <\url{http://statr.me}>
##' @seealso \code{\link[base]{eigen}()}, \code{\link[base]{svd}()},
##' \code{\link[rARPACK]{eigs}()}.
##'
##' @export
##' @rdname svds
##' @keywords array
##' @examples
##' m = 100
##' n = 20
##' k = 5
##' set.seed(111)
##' A = matrix(rnorm(m * n), m)
##' 
##' svds(A, k)
##' svds(t(A), k, nu = 0, nv = 3)
##' 
##' ## Sparse matrices
##' library(Matrix)
##' A[sample(m * n, m * n / 2)] = 0
##' Asp1 = as(A, "dgCMatrix")
##' Asp2 = as(A, "dgRMatrix")
##' 
##' svds(Asp1, k)
##' svds(Asp2, k, nu = 0, nv = 0)
##'
svds <- function(A, k, nu = k, nv = k, opts = list(), ...)
    UseMethod("svds")

##' @rdname svds
##' @export
svds.matrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
{
    fun = if(isSymmetric(A)) svds.real_sym else svds.real_gen
    fun(A, k, nu, nv, opts, ..., mattype = "matrix")
}

##' @rdname svds
##' @export
svds.dgeMatrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
{
    fun = if(isSymmetric(A)) svds.real_sym else svds.real_gen
    fun(A, k, nu, nv, opts, ..., mattype = "dgeMatrix")
}

##' @rdname svds
##' @export
svds.dgCMatrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
{
    fun = if(isSymmetric(A)) svds.real_sym else svds.real_gen
    fun(A, k, nu, nv, opts, ..., mattype = "dgCMatrix")
}

##' @rdname svds
##' @export
svds.dgRMatrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
    svds.real_gen(A, k, nu, nv, opts, ..., mattype = "dgRMatrix")

##' @rdname svds
##' @export
svds.dsyMatrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
    svds.real_sym(A, k, nu, nv, opts, ..., mattype = "dsyMatrix",
                  extra_args = list(use_lower = (A@uplo == "L")))
