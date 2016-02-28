\name{NEWS}
\title{News for Package "rARPACK"}

\section{Changes in rARPACK version 0.10-0}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Updated the backend Spectra library, which fixed the compatibility
            with Eigen >= 3.2.6
    }
  }
}

\section{Changes in rARPACK version 0.9-0}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a bug that causes the algorithm not converging on some matrices
    }
  }
}

\section{Changes in rARPACK version 0.8-1}{
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a compilation problem on Solaris
    }
  }
}

\section{Changes in rARPACK version 0.8-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The backend program is now changed from ARPACK to
            \href{https://github.com/yixuan/spectra}{Spectra},
            which brings cleaner code and better performance.
      \item \code{eigs_sym()} now accepts more matrix types.
      \item Added a C interface for other packages to link to.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed a bug reported by xshi19
            (\href{https://github.com/yixuan/rARPACK/issues/8}{#8})
      \item Fixed a performance issue reported by swajnautcz
            (\href{https://github.com/yixuan/rARPACK/issues/10}{#10})
    }
  }
}

\section{Changes in rARPACK version 0.7-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Support for implicit matrix, contributed by Jiali Mei.
            User can supply a function \code{FUN} rather than an explicit
            matrix to \code{eigs()}, and the eigenvalues/eigenvectors of this
            operator will be computed. \code{FUN(x, args)} must return a vector of the same length as \code{x}.
      \item \code{eigs()} will test the symmetry of matrix before actual
            computation, since symmetric matrices can guarantee real
            eigenvalues and eigenvectors, and the numerical result is more
            stable.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item C++ code of \code{svds()} is completely rewritten. Now it is more
            readable and easier to maintain.
      \item Fix a bug possibly coming from ARPACK, which sometimes gives
            incorrect result of complex eigenvectors.
      \item Avoid using a C random number generator.
    }
  }
}

\section{Changes in rARPACK version 0.6-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Add support for new matrix types: \strong{dgeMatrix} and
            \strong{dgRMatrix}.
      \item \code{eigs()} now allows a full Eigen Decomposition, meaning that
            all the eigenvalues are calculated. In this case \code{eigs()} is
            simply a wrapper of \code{eigen()}, and with a warning issued.
      \item Ditto for \code{svds()}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Rewrite C++ code using classes and templates.
      \item Fix errors in checking the values of \code{k} and \code{ncv}.
    }
  }
}

\section{Changes in rARPACK version 0.5-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Add \code{svds()} function to calculate truncated SVD.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Now sort eigenvalues in decreasing order.
      \item Rename \code{eigs.sym()} to \code{eigs_sym()} to avoid confusion.
      \item Fix a matrix out-of-bound error.
    }
  }
}

\section{Changes in rARPACK version 0.4-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Implement shift-and-invert mode for all supported eigen problems.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Update arpack-ng to 3.1.4.
    }
  }
}

\section{Changes in rARPACK version 0.3-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Now \code{eigs()} supports real symmetric matrices.
    }
  }
}

\section{Changes in rARPACK version 0.2-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Now \code{eigs()} supports sparse real nonsymmetric matrices of the
            class \strong{dgCMatrix}, defined in the \pkg{Matrix} package.
    }
  }
}

\section{Changes in rARPACK version 0.1-0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Initial version. For now \code{eigs()} supports dense real
            nonsymmetric matrices.
    }
  }
}
