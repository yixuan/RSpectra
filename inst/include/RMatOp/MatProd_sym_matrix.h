#ifndef MATPROD_SYM_MATRIX_H
#define MATPROD_SYM_MATRIX_H

#include <Rcpp.h>
#include <R_ext/BLAS.h>  // for BLAS and F77_CALL
#ifndef FCONE
# define FCONE
#endif
#include "MatProd.h"

class MatProd_sym_matrix: public MatProd
{
private:
    const double* mat_ptr;
    const int     n;
    const char    uplo;
    const double  BLAS_alpha;
    const int     BLAS_one;
    const double  BLAS_zero;

public:
    MatProd_sym_matrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        mat_ptr(REAL(mat_)),
        n(nrow_),
        uplo(uplo_),
        BLAS_alpha(1.0),
        BLAS_one(1),
        BLAS_zero(0.0)
    {}

    int rows() const { return n; }
    int cols() const { return n; }

    // y_out = A * x_in
    void perform_op(const double* x_in, double* y_out)
    {
        F77_CALL(dsymv)(&uplo, &n,
                        &BLAS_alpha, mat_ptr, &n,
                        x_in, &BLAS_one, &BLAS_zero,
                        y_out, &BLAS_one FCONE);
    }

    void perform_tprod(const double* x_in, double* y_out)
    {
        perform_op(x_in, y_out);
    }
};


#endif // MATPROD_SYM_MATRIX_H
