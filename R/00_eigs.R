##' Find a Specified Number of Eigenvalues/vectors for Square Matrix
##'
##' @description
##' Given an \eqn{n} by \eqn{n} matrix \eqn{A},
##' function \code{eigs()} can calculate a limited
##' number of eigenvalues and eigenvectors of \eqn{A}.
##' Users can specify the selection criteria by argument
##' \code{which}, e.g., choosing the \eqn{k} largest or smallest
##' eigenvalues and the corresponding eigenvectors.
##' 
##' Currently \code{eigs()} supports matrices of the following classes:
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
##'                          package.\cr
##'   \code{function}   \tab Implicitly specify the matrix through a
##'                          function that has the effect of calculating
##'                          \eqn{f(x)=Ax}{f(x) = A * x}. See section
##'                          \strong{Function Interface} for details.
##' }
##' 
##' \code{eigs_sym()} assumes the matrix is symmetric,
##' and only the lower triangle (or upper triangle, which is
##' controlled by the argument \code{lower}) is used for
##' computation, which guarantees that the eigenvalues and eigenvectors are
##' real, and in some cases reduces the workload. One exception is when
##' \code{A} is a function, in which case the user is responsible for the
##' symmetry of the operator.
##' 
##' \code{eigs_sym()} supports "matrix", "dgeMatrix", "dgCMatrix", "dgRMatrix"
##' and "function" typed matrices.
##' 
##' @param A The matrix whose eigenvalues/vectors are to be computed.
##'          It can also be a function which receives a vector \eqn{x}
##'          and calculates \eqn{Ax}{A * x}.
##'          See section \strong{Function Interface} for details.
##' @param k Number of eigenvalues requested.
##' @param which Selection criteria. See \strong{Details} below.
##' @param sigma Shift parameter. See section \strong{Shift-And-Invert Mode}.
##' @param opts Control parameters related to the computing
##'             algorithm. See \strong{Details} below.
##' @param \dots Currently not used.
##' @param lower For symmetric matrices, should the lower triangle
##'              or upper triangle be used.
##' @param n Only used when \code{A} is a function, to specify the
##'          dimension of the implicit matrix. See section
##'          \strong{Function Interface} for details.
##' @param args Only used when \code{A} is a function. This argument
##'             will be passed to the \code{A} function containing any
##'             extra data. See section \strong{Function Interface}
##'             for details.
##'
##' @details The \code{which} argument is a character string
##' that specifies the type of eigenvalues to be computed.
##' Possible values are:
##'
##' \tabular{ll}{
##'   "LM"  \tab  The \eqn{k} eigenvalues with largest magnitude. Here the
##'               magnitude means the Euclidean norm of complex numbers.\cr
##'   "SM"  \tab  The \eqn{k} eigenvalues with smallest magnitude.\cr
##'   "LR"  \tab  The \eqn{k} eigenvalues with largest real part.\cr
##'   "SR"  \tab  The \eqn{k} eigenvalues with smallest real part.\cr
##'   "LI"  \tab  The \eqn{k} eigenvalues with largest imaginary part.\cr
##'   "SI"  \tab  The \eqn{k} eigenvalues with smallest imaginary part.\cr
##'   "LA"  \tab  The \eqn{k} largest (algebraic) eigenvalues, considering any
##'               negative sign.\cr
##'   "SA"  \tab  The \eqn{k} smallest (algebraic) eigenvalues, considering any
##'               negative sign.\cr
##'   "BE"  \tab  Compute \eqn{k} eigenvalues, half from each end of the
##'               spectrum. When \eqn{k} is odd, compute more from the high
##'               and then from the low end.
##' }
##'
##' \code{eigs()} with matrix type "matrix", "dgeMatrix", "dgCMatrix"
##' and "dgRMatrix" can use "LM",
##' "SM", "LR", "SR", "LI" and "SI".
##' 
##' \code{eigs_sym()}, and \code{eigs()} with matrix type "dsyMatrix"
##' can use "LM", "SM", "LA", "SA" and "BE".
##' 
##' The \code{opts} argument is a list that can supply any of the
##' following parameters:
##'
##' \describe{
##' \item{\code{ncv}}{Number of Lanzcos basis vectors to use. More vectors
##'                   will result in faster convergence, but with greater
##'                   memory use. For general matrix, \code{ncv} must satisfy
##'                   \eqn{k+2\le ncv \le n}{k+2 <= ncv <= n}, and
##'                   for symmetric matrix, the constraint is
##'                   \eqn{k < ncv \le n}{k < ncv <= n}.
##'                   Default is \code{min(n, max(2*k+1, 20))}.}
##' \item{\code{tol}}{Precision parameter. Default is 1e-10.}
##' \item{\code{maxitr}}{Maximum number of iterations. Default is 1000.}
##' \item{\code{retvec}}{Whether to compute eigenvectors. If FALSE,
##'                      only calculate and return eigenvalues.}
##' }
##' 
##' @section Shift-And-Invert Mode:
##' The \code{sigma} argument is used in the shift-and-invert mode.
##' 
##' When \code{sigma} is not \code{NULL}, the selection criteria specified
##' by argument \code{which} will apply to
##' 
##' \deqn{\frac{1}{\lambda-\sigma}}{1/(\lambda-\sigma)}
##' 
##' where \eqn{\lambda}'s are the eigenvalues of \eqn{A}. This mode is useful
##' when user wants to find eigenvalues closest to a given number.
##' For example, if \eqn{\sigma=0}, then \code{which = "LM"} will select the
##' largest values of \eqn{1/|\lambda|}, which turns out to select
##' eigenvalues of \eqn{A} that have the smallest magnitude. The result of
##' using \code{which = "LM", sigma = 0} will be the same as
##' \code{which = "SM"}, but the former one is preferable
##' in that ARPACK is good at finding large
##' eigenvalues rather than small ones. More explanation of the
##' shift-and-invert mode can be found in the SciPy document,
##' \url{http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html}.
##' 
##' @section Function Interface:
##' The matrix \eqn{A} can be specified through a function with
##' the definition
##' 
##' \preformatted{function(x, args)
##' {
##'     ## should return A \%*\% x
##' }}
##' 
##' which receives a vector \code{x} as an argument and returns a vector
##' of the same length. The function should have the effect of calculating
##' \eqn{Ax}{A * x}, and extra arguments can be passed in through the
##' \code{args} parameter. In \code{eigs()}, user should also provide
##' the dimension of the implicit matrix through the argument \code{n}.
##' 
##' @return A list of converged eigenvalues and eigenvectors.
##' \item{values}{Computed eigenvalues.}
##' \item{vectors}{Computed eigenvectors. \code{vectors[, j]} corresponds to \code{values[j]}.}
##' \item{nconv}{Number of converged eigenvalues.}
##' \item{niter}{Number of iterations used in the computation.}
##' \item{nops}{Number of matrix operations used in the computation.}
##' @author Yixuan Qiu \url{http://statr.me}
##' 
##'         Jiali Mei \email{vermouthmjl@@gmail.com}
##' @seealso \code{\link[base]{eigen}()}, \code{\link[base]{svd}()},
##'          \code{\link[rARPACK]{svds}()}
##'
##' @export
##' @rdname eigs
##' @keywords array
##' @examples
##' library(Matrix)
##' n = 20
##' k = 5
##' 
##' ## general matrices have complex eigenvalues
##' set.seed(111)
##' A1 = matrix(rnorm(n^2), n)  ## class "matrix"
##' A2 = Matrix(A1)             ## class "dgeMatrix"
##' 
##' eigs(A1, k)
##' eigs(A2, k, opts = list(retvec = FALSE))  ## eigenvalues only
##' 
##' ## sparse matrices
##' A1[sample(n^2, n^2 / 2)] = 0
##' A3 = as(A1, "dgCMatrix")
##' A4 = as(A1, "dgRMatrix")
##' 
##' eigs(A3, k)
##' eigs(A4, k)
##' 
##' ## function interface
##' f = function(x, args)
##' {
##'     as.numeric(args %*% x)
##' }
##' eigs(f, k, n = n, args = A3)
##' 
##' ## symmetric matrices have real eigenvalues
##' A5 = crossprod(A1)
##' eigs_sym(A5, k)
##' 
##' ## find the smallest (in absolute value) k eigenvalues of A5
##' eigs_sym(A5, k, which = "SM")
##' 
##' ## another way to do this: use the sigma argument
##' eigs_sym(A5, k, sigma = 0)
##' 
##' ## The results should be the same,
##' ## but the latter method is far more stable on large matrices
eigs <- function(A, k, which = "LM", sigma = NULL,
                 opts = list(), ...)
    UseMethod("eigs")

##' @rdname eigs
##' @export
eigs.matrix <- function(A, k, which = "LM", sigma = NULL,
                        opts = list(), ...)
{
    if(isSymmetric(A) &
           which %in% c("LM", "SM", "LR", "SR") &
           (is.null(sigma) || Im(sigma) == 0))
    {
        if(which == "LR")  which = "LA"
        if(which == "SR")  which = "SA"
        eigs.real_sym(A, nrow(A), k, which, sigma, opts, ..., mattype = "sym_matrix",
                      extra_args = list(use_lower = TRUE))
    } else {
        eigs.real_gen(A, nrow(A), k, which, sigma, opts, ..., mattype = "matrix")
    }
}

##' @rdname eigs
##' @export
eigs.dgeMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
{
    if(isSymmetric(A) &
           which %in% c("LM", "SM", "LR", "SR") &
           (is.null(sigma) || Im(sigma) == 0))
    {
        if(which == "LR")  which = "LA"
        if(which == "SR")  which = "SA"
        eigs.real_sym(A, nrow(A), k, which, sigma, opts, ..., mattype = "sym_dgeMatrix",
                      extra_args = list(use_lower = TRUE))
    } else {
        eigs.real_gen(A, nrow(A), k, which, sigma, opts, ..., mattype = "dgeMatrix")
    }
}

##' @rdname eigs
##' @export
eigs.dgCMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
{
    if(isSymmetric(A) &
           which %in% c("LM", "SM", "LR", "SR") &
           (is.null(sigma) || Im(sigma) == 0))
    {
        if(which == "LR")  which = "LA"
        if(which == "SR")  which = "SA"
        eigs.real_sym(A, nrow(A), k, which, sigma, opts, ..., mattype = "sym_dgCMatrix",
                      extra_args = list(use_lower = TRUE))
    } else {
        eigs.real_gen(A, nrow(A), k, which, sigma, opts, ..., mattype = "dgCMatrix")
    }
}

##' @rdname eigs
##' @export
## isSymmetric() does not support dgRMatrix
eigs.dgRMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
    eigs.real_gen(A, nrow(A), k, which, sigma, opts, ..., mattype = "dgRMatrix")

##' @rdname eigs
##' @export
eigs.dsyMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
    eigs.real_sym(A, nrow(A), k, which, sigma, opts, ..., mattype = "dsyMatrix",
                  extra_args = list(use_lower = (A@uplo == "L")))

##' @rdname eigs
##' @export
eigs.function <- function(A, k, which = "LM", sigma = NULL,
                          opts = list(), ...,
                          n = NULL, args = NULL)
    eigs.real_gen(A, as.integer(n), k, which, sigma, opts, ..., mattype = "function",
                  extra_args = list(fun_args = args))



##' @rdname eigs
##' @usage eigs_sym(A, k, which = "LM", sigma = NULL, opts = list(),
##'    lower = TRUE, ...)
##' @export
eigs_sym <- function(A, k, which = "LM", sigma = NULL, opts = list(),
                     lower = TRUE, ...)
    UseMethod("eigs_sym")

eigs_sym.matrix <- function(A, k, which = "LM", sigma = NULL, opts = list(),
                            lower = TRUE, ...)
{
    eigs.real_sym(A, nrow(A), k, which, sigma, opts, ..., mattype = "sym_matrix",
                  extra_args = list(use_lower = as.logical(lower)))
}

eigs_sym.dgeMatrix <- function(A, k, which = "LM", sigma = NULL, opts = list(),
                               lower = TRUE, ...)
{
    eigs.real_sym(A, nrow(A), k, which, sigma, opts, ..., mattype = "sym_dgeMatrix",
                  extra_args = list(use_lower = as.logical(lower)))
}

eigs_sym.dgCMatrix <- function(A, k, which = "LM", sigma = NULL, opts = list(),
                               lower = TRUE, ...)
{
    eigs.real_sym(A, nrow(A), k, which, sigma, opts, ..., mattype = "sym_dgCMatrix",
                  extra_args = list(use_lower = as.logical(lower)))
}

eigs_sym.dgRMatrix <- function(A, k, which = "LM", sigma = NULL, opts = list(),
                               lower = TRUE, ...)
{
    eigs.real_sym(A, nrow(A), k, which, sigma, opts, ..., mattype = "sym_dgRMatrix",
                  extra_args = list(use_lower = as.logical(lower)))
}

##' @rdname eigs
##' @export
eigs_sym.function <- function(A, k, which = "LM", sigma = NULL, opts = list(),
                              lower = TRUE, ..., n = NULL, args = NULL)
{
    eigs.real_sym(A, as.integer(n), k, which, sigma, opts, ..., mattype = "function",
                  extra_args = list(fun_args = args))
}



## Some enumerations

# Matrix types
MAT_TYPE = c("matrix"    = 0L, "sym_matrix"    = 1L,
             "dgeMatrix" = 2L, "sym_dgeMatrix" = 3L,
             "dsyMatrix" = 4L,
             "dgCMatrix" = 5L, "sym_dgCMatrix" = 6L,
             "dgRMatrix" = 7L, "sym_dgRMatrix" = 8L,
             "function"  = 9L)
# Solver types
SOLVER_TYPE = c("regular" = 0L, "real_shift" = 1L, "complex_shift" = 2L)

# Selection rules
EIGS_RULE = c("LM" = 0L, "LR" = 1L, "LI" = 2L, "LA" = 3L, "SM" = 4L,
              "SR" = 5L, "SI" = 6L, "SA" = 7L, "BE" = 8L)