#ifndef REALSHIFT_DGEMATRIX_H
#define REALSHIFT_DGEMATRIX_H

#include <Rcpp.h>
#include <Rdefines.h>  // for R macros
#include "RealShift_matrix.h"

class RealShift_dgeMatrix: public RealShift_matrix
{
public:
    RealShift_dgeMatrix(SEXP mat_, const int nrow_) :
        RealShift_matrix(GET_SLOT(mat_, Rf_install("x")), nrow_)
    {}
};


#endif // REALSHIFT_DGEMATRIX_H
