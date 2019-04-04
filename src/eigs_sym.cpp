#include <RcppEigen.h>
#include <SymEigs.h>
#include <SpectraC.h>
#include "matops.h"
#include "matops_c.h"

using namespace Spectra;
using Rcpp::as;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

/************************ Macros to generate code ************************/

#define EIG_COMMON_CODE                                                        \
if(user_initvec)                                                               \
    eigs.init(initvec);                                                        \
else                                                                           \
    eigs.init();                                                               \
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
Rcpp::RObject run_eigs_sym(
    MatProd* op, int n, int nev, int ncv, int rule,
    int maxitr, double tol, bool retvec,
    bool user_initvec, const double* initvec
)
{
    Rcpp::RObject evals, evecs;
    int nconv = 0, niter = 0, nops = 0;

    EIG_CODE_GENERATOR(REGULAR, MatProd)

    return Rcpp::List::create(
        Rcpp::Named("values")  = evals,
        Rcpp::Named("vectors") = evecs,
        Rcpp::Named("nconv")   = nconv,
        Rcpp::Named("niter")   = niter,
        Rcpp::Named("nops")    = nops
    );
}



RcppExport SEXP eigs_sym(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
)
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

    const double* initvec = NULL;
    bool user_initvec = as<bool>(params_rcpp["user_initvec"]);
    if(user_initvec)
    {
        Rcpp::NumericVector v0 = params_rcpp["initvec"];
        initvec = v0.begin();
    }

    MatProd* op = get_mat_prod_op(A_mat_r, n, n, params_list_r, mattype);
    Rcpp::RObject res = run_eigs_sym(op, n, nev, ncv, rule,
                                     maxitr, tol, retvec,
                                     user_initvec, initvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Regular mode ************************/



/************************ Shift-and-invert mode ************************/
Rcpp::RObject run_eigs_shift_sym(
    RealShift* op, int n, int nev, int ncv, int rule,
    double sigma, int maxitr, double tol, bool retvec,
    bool user_initvec, const double* initvec
)
{
    Rcpp::RObject evals, evecs;
    int nconv = 0, niter = 0, nops = 0;

    EIG_CODE_GENERATOR(REAL_SHIFT, RealShift)

    return Rcpp::List::create(
        Rcpp::Named("values")  = evals,
        Rcpp::Named("vectors") = evecs,
        Rcpp::Named("nconv")   = nconv,
        Rcpp::Named("niter")   = niter,
        Rcpp::Named("nops")    = nops
    );
}

RcppExport SEXP eigs_shift_sym(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
)
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

    const double* initvec = NULL;
    bool user_initvec = as<bool>(params_rcpp["user_initvec"]);
    if(user_initvec)
    {
        Rcpp::NumericVector v0 = params_rcpp["initvec"];
        initvec = v0.begin();
    }

    RealShift* op = get_real_shift_op_sym(A_mat_r, n, params_list_r, mattype);

    Rcpp::RObject res = run_eigs_shift_sym(op, n, nev, ncv, rule,
                                           sigma, maxitr, tol, retvec,
                                           user_initvec, initvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Shift-and-invert mode ************************/



/************************ C interface ************************/
void eigs_sym_c(
    mat_op op, int n, int k,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals, double* evecs, int* info
)
{
    BEGIN_RCPP

    CMatProd cmat_op(op, n, data);
    Rcpp::List res;
    try {
        res = run_eigs_sym((MatProd*) &cmat_op, n, k, opts->ncv, opts->rule,
                           opts->maxitr, opts->tol, opts->retvec != 0,
                           false, NULL);
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
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals, double* evecs, int* info
)
{
    BEGIN_RCPP

    CRealShift cmat_op(op, n, data);
    Rcpp::List res;
    try {
        res = run_eigs_shift_sym((RealShift*) &cmat_op, n, k, opts->ncv, opts->rule,
                                 sigma, opts->maxitr, opts->tol, opts->retvec != 0,
                                 false, NULL);
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
