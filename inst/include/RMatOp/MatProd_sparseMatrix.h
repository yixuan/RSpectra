#ifndef MATPROD_SPARSEMATRIX_H
#define MATPROD_SPARSEMATRIX_H

#include <RcppEigen.h>
#include "MatProd.h"

template <int Storage>
class MatProd_sparseMatrix: public MatProd
{
private:
    typedef Eigen::SparseMatrix<double, Storage> SpMat;
    typedef Eigen::Map<SpMat> MapSpMat;
    typedef Eigen::Map<const Eigen::VectorXd> MapConstVec;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;

    // Map to Eigen sparse matrix
    MapSpMat  mat;
    const int nrow;
    const int ncol;

public:
    MatProd_sparseMatrix(SEXP mat_, const int nrow_, const int ncol_) :
        mat(Rcpp::as<MapSpMat>(mat_)),
        nrow(nrow_),
        ncol(ncol_)
    {}

    int rows() const { return nrow; }
    int cols() const { return ncol; }

    // y_out = A * x_in
    void perform_op(const double* x_in, double* y_out)
    {
        MapConstVec x(x_in, ncol);
        MapVec y(y_out, nrow);
        y.noalias() = mat * x;
    }

    void perform_tprod(const double* x_in, double* y_out)
    {
        MapConstVec x(x_in, nrow);
        MapVec y(y_out, ncol);
        y.noalias() = mat.transpose() * x;
    }
};

// Operations on "dgCMatrix" class, defined in Matrix package
typedef MatProd_sparseMatrix<Eigen::ColMajor> MatProd_dgCMatrix;

// Operations on "dgRMatrix" class, defined in Matrix package
typedef MatProd_sparseMatrix<Eigen::RowMajor> MatProd_dgRMatrix;


#endif // MATPROD_SPARSEMATRIX_H
