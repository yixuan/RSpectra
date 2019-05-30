#include <RcppEigen.h>

typedef Eigen::Map< Eigen::SparseMatrix<double, Eigen::ColMajor> > MapSpMat;
typedef Eigen::Map< Eigen::SparseMatrix<double, Eigen::RowMajor> > MapSpRMat;


RcppExport SEXP is_sym_dgCMatrix(SEXP mat, SEXP tol)
{
    BEGIN_RCPP

    MapSpMat x = Rcpp::as<MapSpMat>(mat);
    const double eps = Rcpp::as<double>(tol);

    // Early return if not a square matrix
    const int n = x.rows();
    if(x.cols() != n)
        return Rcpp::wrap(false);

    for(int j = 0; j < n; j++)
    {
        for(MapSpMat::InnerIterator it(x, j); it; ++it)
        {
            // Only consider i > j
            const int i = it.row();
            if(i <= j)
                continue;

            if(std::abs(it.value() - x.coeff(j, i)) >= eps)
                return Rcpp::wrap(false);
        }
    }

    return Rcpp::wrap(true);

    END_RCPP
}



RcppExport SEXP is_sym_dgRMatrix(SEXP mat, SEXP tol)
{
    BEGIN_RCPP

    MapSpRMat x = Rcpp::as<MapSpRMat>(mat);
    const double eps = Rcpp::as<double>(tol);

    // Early return if not a square matrix
    const int n = x.rows();
    if(x.cols() != n)
        return Rcpp::wrap(false);

    for(int i = 0; i < n; i++)
    {
        for(MapSpRMat::InnerIterator it(x, i); it; ++it)
        {
            // Only consider i < j
            const int j = it.col();
            if(i >= j)
                continue;

            if(std::abs(it.value() - x.coeff(j, i)) >= eps)
                return Rcpp::wrap(false);
        }
    }

    return Rcpp::wrap(true);

    END_RCPP
}
