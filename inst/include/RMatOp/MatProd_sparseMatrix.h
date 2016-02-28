#ifndef MATPROD_SPARSEMATRIX_H
#define MATPROD_SPARSEMATRIX_H

#include <RcppEigen.h>
#include "MatProd.h"

template <int Storage>
class MatProd_sparseMatrix: public MatProd
{
private:
    typedef Eigen::MappedSparseMatrix<double, Storage> MapSpMat;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;

    // Map to Eigen sparse matrix
    MapSpMat mat;
    const int nrow;
    const int ncol;

public:
    MatProd_sparseMatrix(SEXP mat_, const int nrow_, const int ncol_) :
        mat(Rcpp::as<MapSpMat>(mat_)),
        nrow(nrow_),
        ncol(ncol_)
    {}

    int rows() { return nrow; }
    int cols() { return ncol; }

    // y_out = A * x_in
    void perform_op(double *x_in, double *y_out)
    {
        MapVec x(x_in, ncol);
        MapVec y(y_out, nrow);
        y.noalias() = mat * x;
    }

    void perform_tprod(double *x_in, double *y_out)
    {
        MapVec x(x_in, nrow);
        MapVec y(y_out, ncol);
        y.noalias() = mat.transpose() * x;
    }
};

// Operations on "dgCMatrix" class, defined in Matrix package
typedef MatProd_sparseMatrix<Eigen::ColMajor> MatProd_dgCMatrix;

// Operations on "dgRMatrix" class, defined in Matrix package
typedef MatProd_sparseMatrix<Eigen::RowMajor> MatProd_dgRMatrix;


#endif // MATPROD_SPARSEMATRIX_H
