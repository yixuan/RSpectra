#ifndef MATPROD_FUNCTION_H
#define MATPROD_FUNCTION_H

#include <Rcpp.h>
#include "MatProd.h"

class MatProd_function: public MatProd
{
private:
    Rcpp::Function fun;
    int n;
    Rcpp::RObject args;

public:
    MatProd_function(SEXP mat_, const int nrow_, SEXP args_) :
        fun(mat_),
        n(nrow_),
        args(args_)
    {}

    int rows() { return n; }
    int cols() { return n; }

    // y_out = A * x_in
    void perform_op(double *x_in, double *y_out)
    {
        Rcpp::NumericVector x(n);
        std::copy(x_in, x_in + n, x.begin());

        Rcpp::NumericVector y = fun(x, args);
        if(y.length() != n)
            Rcpp::stop("the provided function should return n elements");

        std::copy(y.begin(), y.end(), y_out);
    }

    // y_out = A' * x_in
    void perform_tprod(double *x_in, double *y_out)
    {
        Rcpp::stop("transpose multiplication not implemented for function-typed matrices");
    }

};


#endif // MATPROD_FUNCTION_H
