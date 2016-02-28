#ifndef REALSHIFT_SYM_DGEMATRIX_H
#define REALSHIFT_SYM_DGEMATRIX_H

#include <Rcpp.h>
#include <Rdefines.h>  // for R macros
#include "RealShift_sym_matrix.h"

class RealShift_sym_dgeMatrix: public RealShift_sym_matrix
{
public:
    RealShift_sym_dgeMatrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        RealShift_sym_matrix(GET_SLOT(mat_, Rf_install("x")), nrow_, uplo_)
    {}
};


#endif // REALSHIFT_SYM_DGEMATRIX_H
