#ifndef MATPROD_SYM_SPARSEMATRIX_H
#define MATPROD_SYM_SPARSEMATRIX_H

#include <RcppEigen.h>
#include "MatProd.h"
#include "SparseMatrixMapping.h"

template <int Storage>
class MatProd_sym_sparseMatrix: public MatProd
{
private:
    typedef Eigen::SparseMatrix<double, Storage> SpMat;
    typedef Eigen::Map<SpMat> MapSpMat;
    typedef Eigen::Map<const Eigen::VectorXd> MapConstVec;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;

    // Map to Eigen sparse matrix
    MapSpMat   mat;
    const int  n;
    const char uplo;

public:
    MatProd_sym_sparseMatrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        mat(map_sparse<Storage>(mat_)),
        n(nrow_),
        uplo(uplo_)
    {}

    int rows() const { return n; }
    int cols() const { return n; }

    // y_out = A * x_in
    void perform_op(const double* x_in, double* y_out)
    {
        MapConstVec x(x_in, n);
        MapVec y(y_out, n);

        if(uplo == 'L')
            y.noalias() = mat.template selfadjointView<Eigen::Lower>() * x;
        else
            y.noalias() = mat.template selfadjointView<Eigen::Upper>() * x;
    }

    void perform_tprod(const double* x_in, double* y_out)
    {
        perform_op(x_in, y_out);
    }
};

// Operations on "dgCMatrix" class, defined in Matrix package
typedef MatProd_sym_sparseMatrix<Eigen::ColMajor> MatProd_sym_dgCMatrix;

// Operations on "dgRMatrix" class, defined in Matrix package
typedef MatProd_sym_sparseMatrix<Eigen::RowMajor> MatProd_sym_dgRMatrix;


#endif // MATPROD_SYM_SPARSEMATRIX_H
