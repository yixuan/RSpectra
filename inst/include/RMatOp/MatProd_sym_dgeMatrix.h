#ifndef MATPROD_SYM_DGEMATRIX_H
#define MATPROD_SYM_DGEMATRIX_H

#include <Rcpp.h>
#include <Rdefines.h>  // for R macros
#include "MatProd_sym_matrix.h"

class MatProd_sym_dgeMatrix: public MatProd_sym_matrix
{
public:
    MatProd_sym_dgeMatrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        MatProd_sym_matrix(GET_SLOT(mat_, Rf_install("x")), nrow_, uplo_)
    {}
};


#endif // MATPROD_SYM_DGEMATRIX_H
