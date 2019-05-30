#ifndef MATPROD_SPARSE_MATRIX_MAPPING_H
#define MATPROD_SPARSE_MATRIX_MAPPING_H

#include <RcppEigen.h>

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


#endif // MATPROD_SPARSE_MATRIX_MAPPING_H
