#ifndef COMPLEXSHIFT_DGEMATRIX_H
#define COMPLEXSHIFT_DGEMATRIX_H

#include <Rcpp.h>
#include <Rdefines.h>  // for R macros
#include "ComplexShift_matrix.h"

class ComplexShift_dgeMatrix: public ComplexShift_matrix
{
public:
    ComplexShift_dgeMatrix(SEXP mat_, const int nrow_) :
        ComplexShift_matrix(GET_SLOT(mat_, Rf_install("x")), nrow_)
    {}
};


#endif // COMPLEXSHIFT_DGEMATRIX_H
