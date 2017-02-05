#ifndef MATPROD_FUNCTION_H
#define MATPROD_FUNCTION_H

#include <Rcpp.h>
#include "MatProd.h"

class MatProd_function: public MatProd
{
private:
    Rcpp::Function A;
    Rcpp::Function Atrans;
    const int      nrow;
    const int      ncol;
    Rcpp::RObject  args;

public:
    MatProd_function(SEXP mat_, SEXP trans_, const int nrow_, const int ncol_, SEXP args_) :
        A(mat_),
        Atrans(trans_),
        nrow(nrow_),
        ncol(ncol_),
        args(args_)
    {}

    int rows() const { return nrow; }
    int cols() const { return ncol; }

    // y_out = A * x_in
    void perform_op(const double* x_in, double* y_out)
    {
        Rcpp::NumericVector x(ncol);
        std::copy(x_in, x_in + ncol, x.begin());

        Rcpp::NumericVector y = A(x, args);
        if(y.length() != nrow)
            Rcpp::stop("the provided function should return m elements");

        std::copy(y.begin(), y.end(), y_out);
    }

    // y_out = A' * x_in
    void perform_tprod(const double* x_in, double* y_out)
    {
        Rcpp::NumericVector x(nrow);
        std::copy(x_in, x_in + nrow, x.begin());

        Rcpp::NumericVector y = Atrans(x, args);
        if(y.length() != ncol)
            Rcpp::stop("the provided transpose function should return n elements");

        std::copy(y.begin(), y.end(), y_out);
    }

};


#endif // MATPROD_FUNCTION_H
