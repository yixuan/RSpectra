#ifndef MATPROD_DGEMATRIX_H
#define MATPROD_DGEMATRIX_H

#include <Rcpp.h>
#include <Rdefines.h>  // for R macros
#include "MatProd_matrix.h"

class MatProd_dgeMatrix: public MatProd_matrix
{
public:
    MatProd_dgeMatrix(SEXP mat_, const int nrow_, const int ncol_) :
        MatProd_matrix(GET_SLOT(mat_, Rf_install("x")), nrow_, ncol_)
    {}
};


#endif // MATPROD_DGEMATRIX_H
