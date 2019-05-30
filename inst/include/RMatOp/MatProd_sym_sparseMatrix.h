#ifndef MATPROD_SYM_SPARSEMATRIX_H
#define MATPROD_SYM_SPARSEMATRIX_H

#include <RcppEigen.h>
#include "MatProd.h"

// Mapping a dgCMatrix/dsCMatrix/dgRMatrix/dsRMatrix R object to Eigen
// Default is ColMajor
template <int Storage>
inline Eigen::Map< Eigen::SparseMatrix<double, Storage> > map_sparse(SEXP mat)
{
    Rcpp::S4 obj(mat);
    if(!(obj.is("dgCMatrix") || obj.is("dsCMatrix")))
        throw std::invalid_argument("Need S4 class dgCMatrix or dsCMatrix for a mapped sparse matrix");

    Rcpp::IntegerVector dim(obj.slot("Dim")), i(obj.slot("i")), p(obj.slot("p"));
    Rcpp::NumericVector x(obj.slot("x"));

    return Eigen::Map< Eigen::SparseMatrix<double, Storage> >(
        dim[0], dim[1], p[dim[1]], p.begin(), i.begin(), x.begin()
    );
}

// Specialization for RowMajor
template <>
inline Eigen::Map< Eigen::SparseMatrix<double, Eigen::RowMajor> > map_sparse<Eigen::RowMajor>(SEXP mat)
{
    Rcpp::S4 obj(mat);
    if(!(obj.is("dgRMatrix") || obj.is("dsRMatrix")))
        throw std::invalid_argument("Need S4 class dgRMatrix or dsRMatrix for a mapped sparse matrix");

    Rcpp::IntegerVector dim(obj.slot("Dim")), j(obj.slot("j")), p(obj.slot("p"));
    Rcpp::NumericVector x(obj.slot("x"));

    return Eigen::Map< Eigen::SparseMatrix<double, Eigen::RowMajor> >(
            dim[0], dim[1], p[dim[1]], p.begin(), j.begin(), x.begin()
    );
}



template <int Storage>
class MatProd_sym_sparseMatrix: public MatProd
{
private:
    typedef Eigen::Map< Eigen::SparseMatrix<double, Storage> > MapSpMat;
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
