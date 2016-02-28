#include <RcppEigen.h>
#include <SymEigs.h>
#include "matops.h"

using namespace Spectra;
using Rcpp::as;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

enum SOLVER_TYPE {
    REGULAR = 0,
    REAL_SHIFT
};

/************************ Macros to generate code ************************/

#define EIG_COMMON_CODE                                                        \
eigs.init(init_resid);                                                         \
nconv = eigs.compute(maxitr, tol);                                             \
if(nconv < nev)                                                                \
    Rcpp::warning("only %d eigenvalue(s) converged, less than k = %d",         \
                  nconv, nev);                                                 \
evals = Rcpp::wrap(eigs.eigenvalues());                                        \
if(retvec)                                                                     \
{                                                                              \
    Rcpp::NumericMatrix evecs_ret(n, nconv);                                   \
    MapMat m(evecs_ret.begin(), n, nconv);                                     \
    m.noalias() = eigs.eigenvectors();                                         \
    evecs = evecs_ret;                                                         \
} else {                                                                       \
    evecs = R_NilValue;                                                        \
}                                                                              \
niter = eigs.num_iterations();                                                 \
nops = eigs.num_operations();



#define EIG_CODE_REGULAR(RULE, OPTYPE)                                         \
SymEigsSolver<double, RULE, OPTYPE> eigs(op, nev, ncv);                        \
EIG_COMMON_CODE



#define EIG_CODE_REAL_SHIFT(RULE, OPTYPE)                                      \
SymEigsShiftSolver<double, RULE, OPTYPE> eigs(op, nev, ncv, sigma);            \
EIG_COMMON_CODE



#define EIG_CODE_GENERATOR(SOLVER, OPTYPE)                                     \
switch(rule)                                                                   \
{                                                                              \
    case WHICH_LM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LM, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_LA :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LA, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_SM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SM, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_SA :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SA, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_BE :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_BE, OPTYPE) }                              \
        break;                                                                 \
    default:                                                                   \
        Rcpp::stop("unsupported selection rule");                              \
}

/************************ Macros to generate code ************************/



/************************ Regular mode ************************/
Rcpp::RObject run_eigs_sym(MatProd* op, int n, int nev, int ncv, int rule,
                           int maxitr, double tol, bool retvec)
{
    // Prepare initial residuals
    double *init_resid;
    #include "rands.h"
    if(n <= rands_len)
    {
        init_resid = rands;
    } else {
        init_resid = new double[n];
        double *coef_pntr = init_resid;
        for(int i = 0; i < n / rands_len; i++, coef_pntr += rands_len)
        {
            std::copy(rands, rands + rands_len, coef_pntr);
        }
        std::copy(rands, rands + n % rands_len, coef_pntr);
    }

    Rcpp::RObject evals, evecs;
    int nconv = 0, niter = 0, nops = 0;

    EIG_CODE_GENERATOR(REGULAR, MatProd)

    if(n > rands_len)
        delete [] init_resid;

    return Rcpp::List::create(
        Rcpp::Named("values")  = evals,
        Rcpp::Named("vectors") = evecs,
        Rcpp::Named("nconv")   = nconv,
        Rcpp::Named("niter")   = niter,
        Rcpp::Named("nops")    = nops
    );
}



RcppExport SEXP eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);

    int n        = as<int>(n_scalar_r);
    int nev      = as<int>(k_scalar_r);
    int ncv      = as<int>(params_rcpp["ncv"]);
    int rule     = as<int>(params_rcpp["which"]);
    double tol   = as<double>(params_rcpp["tol"]);
    int maxitr   = as<int>(params_rcpp["maxitr"]);
    bool retvec  = as<bool>(params_rcpp["retvec"]);
    int mattype  = as<int>(mattype_scalar_r);

    MatProd *op = get_mat_prod_op(A_mat_r, n, n, params_list_r, mattype);
    Rcpp::RObject res = run_eigs_sym(op, n, nev, ncv, rule, maxitr, tol, retvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Regular mode ************************/



/************************ Shift-and-invert mode ************************/
Rcpp::RObject run_eigs_shift_sym(RealShift* op, int n, int nev, int ncv, int rule,
                                 double sigma, int maxitr, double tol, bool retvec)
{
    // Prepare initial residuals
    double *init_resid;
    #include "rands.h"
    if(n <= rands_len)
    {
        init_resid = rands;
    } else {
        init_resid = new double[n];
        double *coef_pntr = init_resid;
        for(int i = 0; i < n / rands_len; i++, coef_pntr += rands_len)
        {
            std::copy(rands, rands + rands_len, coef_pntr);
        }
        std::copy(rands, rands + n % rands_len, coef_pntr);
    }

    Rcpp::RObject evals, evecs;
    int nconv = 0, niter = 0, nops = 0;

    EIG_CODE_GENERATOR(REAL_SHIFT, RealShift)

    if(n > rands_len)
        delete [] init_resid;

    return Rcpp::List::create(
        Rcpp::Named("values")  = evals,
        Rcpp::Named("vectors") = evecs,
        Rcpp::Named("nconv")   = nconv,
        Rcpp::Named("niter")   = niter,
        Rcpp::Named("nops")    = nops
    );
}

RcppExport SEXP eigs_shift_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                               SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);

    int n        = as<int>(n_scalar_r);
    int nev      = as<int>(k_scalar_r);
    int ncv      = as<int>(params_rcpp["ncv"]);
    int rule     = as<int>(params_rcpp["which"]);
    double tol   = as<double>(params_rcpp["tol"]);
    int maxitr   = as<int>(params_rcpp["maxitr"]);
    bool retvec  = as<bool>(params_rcpp["retvec"]);
    int mattype  = as<int>(mattype_scalar_r);
    double sigma = as<double>(params_rcpp["sigma"]);

    RealShift *op = eigs_sym_get_real_shift_op(A_mat_r, n, params_list_r, mattype);

    Rcpp::RObject res = run_eigs_shift_sym(op, n, nev, ncv, rule, sigma, maxitr, tol, retvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Shift-and-invert mode ************************/



/************************ C interface ************************/
#include "c_interface.h"

void eigs_sym_c(
    mat_op op, int n, int k,
    const arpack_opts *opts, void *data,
    int *nconv, int *niter, int *nops,
    double *evals, double *evecs, int *info
)
{
    BEGIN_RCPP

    CMatProd cmat_op(op, n, data);
    Rcpp::List res;
    try {
        res = run_eigs_sym((MatProd*) &cmat_op, n, k, opts->ncv, opts->rule,
                           opts->maxitr, opts->tol, opts->retvec != 0);
        *info = 0;
    } catch(...) {
        *info = 1;  // indicates error
    }

    *nconv = Rcpp::as<int>(res["nconv"]);
    *niter = Rcpp::as<int>(res["niter"]);
    *nops  = Rcpp::as<int>(res["nops"]);
    Rcpp::NumericVector val = res["values"];
    std::copy(val.begin(), val.end(), evals);
    if(opts->retvec != 0)
    {
        Rcpp::NumericMatrix vec = res["vectors"];
        std::copy(vec.begin(), vec.end(), evecs);
    }

    VOID_END_RCPP
}

void eigs_sym_shift_c(
    mat_op op, int n, int k, double sigma,
    const arpack_opts *opts, void *data,
    int *nconv, int *niter, int *nops,
    double *evals, double *evecs, int *info
)
{
    BEGIN_RCPP

    CRealShift cmat_op(op, n, data);
    Rcpp::List res;
    try {
        res = run_eigs_shift_sym((RealShift*) &cmat_op, n, k, opts->ncv, opts->rule,
                                 sigma, opts->maxitr, opts->tol, opts->retvec != 0);
        *info = 0;
    } catch(...) {
        *info = 1;  // indicates error
    }

    *nconv = Rcpp::as<int>(res["nconv"]);
    *niter = Rcpp::as<int>(res["niter"]);
    *nops  = Rcpp::as<int>(res["nops"]);
    Rcpp::NumericVector val = res["values"];
    std::copy(val.begin(), val.end(), evals);
    if(opts->retvec != 0)
    {
        Rcpp::NumericMatrix vec = res["vectors"];
        std::copy(vec.begin(), vec.end(), evecs);
    }

    VOID_END_RCPP
}
/************************ C interface ************************/
