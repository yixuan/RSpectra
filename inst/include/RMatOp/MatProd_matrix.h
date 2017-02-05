#ifndef MATPROD_MATRIX_H
#define MATPROD_MATRIX_H

// #include <RcppEigen.h>
#include <R_ext/BLAS.h>  // for BLAS and F77_CALL
#include "MatProd.h"

class MatProd_matrix: public MatProd
{
private:
//  typedef Eigen::Map<const Eigen::MatrixXd> MapMat;
//  typedef Eigen::Map<Eigen::VectorXd> MapVec;

    const double* mat_ptr;
//  const MapMat  mat;
    const int     nrow;
    const int     ncol;

    const double  BLAS_alpha;
    const int     BLAS_one;
    const double  BLAS_zero;


public:
    MatProd_matrix(SEXP mat_, const int nrow_, const int ncol_) :
        mat_ptr(REAL(mat_)),
//      mat(mat_ptr, nrow_, ncol_),
        nrow(nrow_),
        ncol(ncol_),

        BLAS_alpha(1.0),
        BLAS_one(1),
        BLAS_zero(0.0)

    {}

    int rows() const { return nrow; }
    int cols() const { return ncol; }

    // y_out = A * x_in
    void perform_op(const double* x_in, double* y_out)
    {
/*
        MapVec x(x_in, ncol);
        MapVec y(y_out, nrow);
        y.noalias() = mat * x;
*/
        F77_CALL(dgemv)("N", &nrow, &ncol,
                        &BLAS_alpha, mat_ptr, &nrow,
                        x_in, &BLAS_one, &BLAS_zero,
                        y_out, &BLAS_one);
    }

    void perform_tprod(const double* x_in, double* y_out)
    {
/*
        MapVec x(x_in, nrow);
        MapVec y(y_out, ncol);
        y.noalias() = mat.transpose() * x;
*/
        F77_CALL(dgemv)("T", &nrow, &ncol,
                        &BLAS_alpha, mat_ptr, &nrow,
                        x_in, &BLAS_one, &BLAS_zero,
                        y_out, &BLAS_one);
    }
};


#endif // MATPROD_MATRIX_H
