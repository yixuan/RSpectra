#include <RcppEigen.h>
#include "matops.h"

MatProd* get_mat_prod_op(SEXP mat, int nrow, int ncol, SEXP extra_arg, int mat_type)
{
    MatProd* op;

    Rcpp::List args(extra_arg);

    switch(mat_type)
    {
    case MATRIX:
        op = new MatProd_matrix(mat, nrow, ncol);
        break;
    case SYM_MATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new MatProd_sym_matrix(mat, nrow, use_lower ? 'L' : 'U');
        }
        break;
    case DGEMATRIX:
        op = new MatProd_dgeMatrix(mat, nrow, ncol);
        break;
    case SYM_DGEMATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new MatProd_sym_dgeMatrix(mat, nrow, use_lower ? 'L' : 'U');
        }
        break;
    case DSYMATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new MatProd_dsyMatrix(mat, nrow, use_lower ? 'L' : 'U');
        }
        break;
    case DGCMATRIX:
        op = new MatProd_dgCMatrix(mat, nrow, ncol);
        break;
    case SYM_DGCMATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new MatProd_sym_dgCMatrix(mat, nrow, use_lower ? 'L' : 'U');
        }
        break;
    case DGRMATRIX:
        op = new MatProd_dgRMatrix(mat, nrow, ncol);
        break;
    case SYM_DGRMATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new MatProd_sym_dgRMatrix(mat, nrow, use_lower ? 'L' : 'U');
        }
        break;
    case FUNCTION:
        {
        SEXP trans    = args["Atrans"];
        SEXP fun_args = args["fun_args"];
        op = new MatProd_function(mat, trans, nrow, ncol, fun_args);
        }
        break;
    default:
        Rcpp::stop("unsupported matrix type");
        // Eliminate compiler warning, but should not reach here
        op = new MatProd_matrix(mat, nrow, ncol);
    }

    return op;
}

RealShift* get_real_shift_op_sym(SEXP mat, int n, SEXP extra_arg, int mat_type)
{
    RealShift* op;

    Rcpp::List args(extra_arg);

    switch(mat_type)
    {
    case MATRIX:
        op = new RealShift_matrix(mat, n);
        break;
    case SYM_MATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new RealShift_sym_matrix(mat, n, use_lower ? 'L' : 'U');
        }
        break;
    case DGEMATRIX:
        op = new RealShift_dgeMatrix(mat, n);
        break;
    case SYM_DGEMATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new RealShift_sym_dgeMatrix(mat, n, use_lower ? 'L' : 'U');
        }
        break;
    case DSYMATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new RealShift_dsyMatrix(mat, n, use_lower ? 'L' : 'U');
        }
        break;
    case DGCMATRIX:
        op = new RealShift_dgCMatrix(mat, n);
        break;
    case SYM_DGCMATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new RealShift_sym_dgCMatrix(mat, n, use_lower ? 'L' : 'U');
        }
        break;
    case DGRMATRIX:
        op = new RealShift_dgRMatrix(mat, n);
        break;
    case SYM_DGRMATRIX:
        {
        bool use_lower = Rcpp::as<bool>(args["use_lower"]);
        op = new RealShift_sym_dgRMatrix(mat, n, use_lower ? 'L' : 'U');
        }
        break;
    default:
        Rcpp::stop("unsupported matrix type");
        // Eliminate compiler warning, but should not reach here
        op = new RealShift_matrix(mat, n);
    }

    return op;
}

RealShift* get_real_shift_op_gen(SEXP mat, int n, SEXP extra_arg, int mat_type)
{
    RealShift* op;

    Rcpp::List args(extra_arg);

    switch(mat_type)
    {
    case MATRIX:
        op = new RealShift_matrix(mat, n);
        break;
    case DGEMATRIX:
        op = new RealShift_dgeMatrix(mat, n);
        break;
    case DGCMATRIX:
        op = new RealShift_dgCMatrix(mat, n);
        break;
    case DGRMATRIX:
        op = new RealShift_dgRMatrix(mat, n);
        break;
    default:
        Rcpp::stop("unsupported matrix type");
        // Eliminate compiler warning, but should not reach here
        op = new RealShift_matrix(mat, n);
    }

    return op;
}

ComplexShift* get_complex_shift_op(SEXP mat, int n, SEXP extra_arg, int mat_type)
{
    ComplexShift* op;

    Rcpp::List args(extra_arg);

    switch(mat_type)
    {
    case MATRIX:
        op = new ComplexShift_matrix(mat, n);
        break;
    case DGEMATRIX:
        op = new ComplexShift_dgeMatrix(mat, n);
        break;
    case DGCMATRIX:
        op = new ComplexShift_dgCMatrix(mat, n);
        break;
    case DGRMATRIX:
        op = new ComplexShift_dgRMatrix(mat, n);
        break;
    default:
        Rcpp::stop("unsupported matrix type");
        // Eliminate compiler warning, but should not reach here
        op = new ComplexShift_matrix(mat, n);
    }

    return op;
}
