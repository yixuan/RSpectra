#ifndef MATOPS_H
#define MATOPS_H

#include <RMatOp.h>
#include <Rinternals.h>

enum MAT_TYPE {
    MATRIX = 0,
    SYM_MATRIX,
    DGEMATRIX,
    SYM_DGEMATRIX,
    DSYMATRIX,
    DGCMATRIX,
    SYM_DGCMATRIX,
    DGRMATRIX,
    SYM_DGRMATRIX,
    FUNCTION
};

MatProd* get_mat_prod_op(SEXP mat, int nrow, int ncol, SEXP extra_arg, int mat_type);

RealShift* eigs_sym_get_real_shift_op(SEXP mat, int n, SEXP extra_arg, int mat_type);

RealShift* eigs_gen_get_real_shift_op(SEXP mat, int n, SEXP extra_arg, int mat_type);

ComplexShift* get_complex_shift_op(SEXP mat, int n, SEXP extra_arg, int mat_type);


#endif // MATOPS_H
