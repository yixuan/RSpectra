#include <RcppEigen.h>
#include <GenEigs.h>
#include <SpectraC.h>
#include "matops.h"
#include "matops_c.h"

using namespace Spectra;
using Rcpp::as;

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
    evecs = Rcpp::wrap(eigs.eigenvectors());                                   \
else                                                                           \
    evecs = R_NilValue;                                                        \
niter = eigs.num_iterations();                                                 \
nops = eigs.num_operations();



#define EIG_CODE_REGULAR(RULE, OPTYPE)                                         \
GenEigsSolver<double, RULE, OPTYPE> eigs(op, nev, ncv);                        \
EIG_COMMON_CODE



#define EIG_CODE_REAL_SHIFT(RULE, OPTYPE)                                      \
GenEigsRealShiftSolver<double, RULE, OPTYPE> eigs(op, nev, ncv, sigmar);       \
EIG_COMMON_CODE



#define EIG_CODE_COMPLEX_SHIFT(RULE, OPTYPE)                                   \
GenEigsComplexShiftSolver<double, RULE, OPTYPE> eigs(op, nev, ncv, sigmar, sigmai);      \
EIG_COMMON_CODE



#define EIG_CODE_GENERATOR(SOLVER, OPTYPE)                                     \
switch(rule)                                                                   \
{                                                                              \
    case WHICH_LM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LM, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_LR :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LR, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_LI :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LI, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_SM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SM, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_SR :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SR, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_SI :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SI, OPTYPE) }                              \
        break;                                                                 \
    default:                                                                   \
        Rcpp::stop("unsupported selection rule");                              \
}

/************************ Macros to generate code ************************/



/************************ Regular mode ************************/
Rcpp::RObject run_eigs_gen(
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


RcppExport SEXP eigs_gen(
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

    MatProd* op = get_mat_prod_op(A_mat_r, n, n, params_list_r, mattype);;
    Rcpp::RObject res = run_eigs_gen(op, n, nev, ncv, rule,
                                     maxitr, tol, retvec,
                                     user_initvec, initvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Regular mode ************************/



/************************ Real shift mode ************************/
Rcpp::RObject run_eigs_real_shift_gen(
    RealShift* op, int n, int nev, int ncv, int rule,
    double sigmar, int maxitr, double tol, bool retvec,
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

RcppExport SEXP eigs_real_shift_gen(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);

    int n         = as<int>(n_scalar_r);
    int nev       = as<int>(k_scalar_r);
    int ncv       = as<int>(params_rcpp["ncv"]);
    int rule      = as<int>(params_rcpp["which"]);
    double tol    = as<double>(params_rcpp["tol"]);
    int maxitr    = as<int>(params_rcpp["maxitr"]);
    bool retvec   = as<bool>(params_rcpp["retvec"]);
    int mattype   = as<int>(mattype_scalar_r);
    double sigmar = as<double>(params_rcpp["sigmar"]);

    const double* initvec = NULL;
    bool user_initvec = as<bool>(params_rcpp["user_initvec"]);
    if(user_initvec)
    {
        Rcpp::NumericVector v0 = params_rcpp["initvec"];
        initvec = v0.begin();
    }

    RealShift* op = get_real_shift_op_gen(A_mat_r, n, params_list_r, mattype);
    Rcpp::RObject res = run_eigs_real_shift_gen(op, n, nev, ncv, rule,
                                                sigmar, maxitr, tol, retvec,
                                                user_initvec, initvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Real shift mode ************************/



/************************ Complex shift mode ************************/
Rcpp::RObject run_eigs_complex_shift_gen(
    ComplexShift* op, int n, int nev, int ncv, int rule,
    double sigmar, double sigmai, int maxitr, double tol, bool retvec,
    bool user_initvec, const double* initvec
)
{
    Rcpp::RObject evals, evecs;
    int nconv = 0, niter = 0, nops = 0;

    EIG_CODE_GENERATOR(COMPLEX_SHIFT, ComplexShift)

    return Rcpp::List::create(
        Rcpp::Named("values")  = evals,
        Rcpp::Named("vectors") = evecs,
        Rcpp::Named("nconv")   = nconv,
        Rcpp::Named("niter")   = niter,
        Rcpp::Named("nops")    = nops
    );
}

RcppExport SEXP eigs_complex_shift_gen(
    SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
    SEXP params_list_r, SEXP mattype_scalar_r
)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);

    int n         = as<int>(n_scalar_r);
    int nev       = as<int>(k_scalar_r);
    int ncv       = as<int>(params_rcpp["ncv"]);
    int rule      = as<int>(params_rcpp["which"]);
    double tol    = as<double>(params_rcpp["tol"]);
    int maxitr    = as<int>(params_rcpp["maxitr"]);
    bool retvec   = as<bool>(params_rcpp["retvec"]);
    int mattype   = as<int>(mattype_scalar_r);
    double sigmar = as<double>(params_rcpp["sigmar"]);
    double sigmai = as<double>(params_rcpp["sigmai"]);

    const double* initvec = NULL;
    bool user_initvec = as<bool>(params_rcpp["user_initvec"]);
    if(user_initvec)
    {
        Rcpp::NumericVector v0 = params_rcpp["initvec"];
        initvec = v0.begin();
    }

    ComplexShift* op = get_complex_shift_op(A_mat_r, n, params_list_r, mattype);
    Rcpp::RObject res = run_eigs_complex_shift_gen(op, n, nev, ncv, rule, sigmar, sigmai,
                                                   maxitr, tol, retvec,
                                                   user_initvec, initvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Complex shift mode ************************/



/************************ C interface ************************/
void eigs_gen_c(
    mat_op op, int n, int k,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals_r, double* evals_i, double* evecs_r, double* evecs_i, int* info
)
{
    BEGIN_RCPP

    CMatProd cmat_op(op, n, data);
    Rcpp::List res;
    try {
        res = run_eigs_gen((MatProd*) &cmat_op, n, k, opts->ncv, opts->rule,
                           opts->maxitr, opts->tol, opts->retvec != 0,
                           false, NULL);
        *info = 0;
    } catch(...) {
        *info = 1;  // indicates error
    }

    *nconv = Rcpp::as<int>(res["nconv"]);
    *niter = Rcpp::as<int>(res["niter"]);
    *nops  = Rcpp::as<int>(res["nops"]);
    Rcpp::ComplexVector val = res["values"];
    for(int i = 0; i < val.length(); i++)
    {
        Rcomplex c = val[i];
        evals_r[i] = c.r;
        evals_i[i] = c.i;
    }
    if(opts->retvec != 0)
    {
        Rcpp::ComplexMatrix vec = res["vectors"];
        for(int i = 0; i < vec.length(); i++)
        {
            Rcomplex c = val[i];
            evecs_r[i] = c.r;
            evecs_i[i] = c.i;
        }
    }

    VOID_END_RCPP
}

void eigs_gen_shift_c(
    mat_op op, int n, int k, double sigmar, double sigmai,
    const spectra_opts* opts, void* data,
    int* nconv, int* niter, int* nops,
    double* evals_r, double* evals_i, double* evecs_r, double* evecs_i, int* info
)
{
    BEGIN_RCPP

    Rcpp::List res;
    try {
        if(std::abs(sigmai) <= 1e-12)
        {
            CRealShift cmat_op(op, n, data);
            res = run_eigs_real_shift_gen(
                      (RealShift*) &cmat_op, n, k, opts->ncv, opts->rule,
                      sigmar, opts->maxitr, opts->tol, opts->retvec != 0,
                      false, NULL
                  );
        }
        else {
            CComplexShift cmat_op(op, n, data);
            res = run_eigs_complex_shift_gen(
                      (ComplexShift*) &cmat_op, n, k, opts->ncv, opts->rule,
                      sigmar, sigmai, opts->maxitr, opts->tol, opts->retvec != 0,
                      false, NULL
                  );
        }

        *info = 0;
    } catch(...) {
        *info = 1;  // indicates error
    }

    *nconv = Rcpp::as<int>(res["nconv"]);
    *niter = Rcpp::as<int>(res["niter"]);
    *nops  = Rcpp::as<int>(res["nops"]);
    Rcpp::ComplexVector val = res["values"];
    for(int i = 0; i < val.length(); i++)
    {
        Rcomplex c = val[i];
        evals_r[i] = c.r;
        evals_i[i] = c.i;
    }
    if(opts->retvec != 0)
    {
        Rcpp::ComplexMatrix vec = res["vectors"];
        for(int i = 0; i < vec.length(); i++)
        {
            Rcomplex c = val[i];
            evecs_r[i] = c.r;
            evecs_i[i] = c.i;
        }
    }

    VOID_END_RCPP
}
/************************ C interface ************************/
