#ifndef REALSHIFT_DSYMATRIX_H
#define REALSHIFT_DSYMATRIX_H

#include <Rcpp.h>
#include <Rdefines.h>  // for R macros
#include "RealShift_sym_matrix.h"

class RealShift_dsyMatrix : public RealShift_sym_matrix
{
public:
    RealShift_dsyMatrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        RealShift_sym_matrix(GET_SLOT(mat_, Rf_install("x")), nrow_, uplo_)
    {}
};


#endif // REALSHIFT_DSYMATRIX_H
