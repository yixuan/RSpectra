#include <RcppEigen.h>
#include <SymEigs.h>
#include "matops.h"

using namespace Spectra;
using Rcpp::as;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

RcppExport SEXP svds_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP nu_scalar_r, SEXP nv_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_svds(params_list_r);

    int n        = as<int>(n_scalar_r);
    int k        = as<int>(k_scalar_r);
    int nu       = as<int>(nu_scalar_r);
    int nv       = as<int>(nv_scalar_r);
    int ncv      = as<int>(params_svds["ncv"]);
    double tol   = as<double>(params_svds["tol"]);
    int maxitr   = as<int>(params_svds["maxitr"]);
    int mattype  = as<int>(mattype_scalar_r);

    MatProd *op = get_mat_prod_op(A_mat_r, n, n, params_list_r, mattype);

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

    SymEigsSolver<double, LARGEST_MAGN, MatProd> eigs(op, k, ncv);
    eigs.init(init_resid);
    int nconv = eigs.compute(maxitr, tol);
    if(nconv < k)
        Rcpp::warning("only %d singular values converged, less than k = %d", nconv, k);

    nu = std::min(nu, nconv);
    nv = std::min(nv, nconv);

    Eigen::VectorXd evals = eigs.eigenvalues();
    Eigen::MatrixXd evecs = eigs.eigenvectors();

    Rcpp::NumericVector d(nconv);
    Rcpp::NumericMatrix u(n, nu), v(n, nv);

    // Copy evals to d
    std::copy(evals.data(), evals.data() + nconv, d.begin());

    // Copy evecs to u and v with the specified number of columns
    std::copy(evecs.data(), evecs.data() + nu * n, u.begin());
    std::copy(evecs.data(), evecs.data() + nv * n, v.begin());

    // We need to make sure that singular values are nonnegative,
    // so move the sign to v.
    for(int i = 0; i < nconv; i++)
    {
        if(d[i] < 0)
        {
            d[i] = -d[i];
            if(i < nv)
            {
                double *ptr = &v(0, i);
                std::transform(ptr, ptr + n, ptr, std::negate<double>());
            }
        }
    }

    Rcpp::RObject u_ret;
    Rcpp::RObject v_ret;

    if(nu > 0)  u_ret = u;  else  u_ret = R_NilValue;
    if(nv > 0)  v_ret = v;  else  v_ret = R_NilValue;

    if(n > rands_len)
        delete [] init_resid;

    delete op;

    return Rcpp::List::create(
        Rcpp::Named("d")     = d,
        Rcpp::Named("u")     = u_ret,
        Rcpp::Named("v")     = v_ret,
        Rcpp::Named("niter") = eigs.num_iterations(),
        Rcpp::Named("nops")  = eigs.num_operations()
    );

    END_RCPP
}



inline double simple_sqrt(double x) { return std::sqrt(x); }

RcppExport SEXP svds_gen(SEXP A_mat_r, SEXP m_scalar_r, SEXP n_scalar_r,
                         SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_svds(params_list_r);

    int m        = as<int>(m_scalar_r);
    int n        = as<int>(n_scalar_r);
    int k        = as<int>(k_scalar_r);
    int nu       = as<int>(nu_scalar_r);
    int nv       = as<int>(nv_scalar_r);
    int ncv      = as<int>(params_svds["ncv"]);
    double tol   = as<double>(params_svds["tol"]);
    int maxitr   = as<int>(params_svds["maxitr"]);
    int mattype  = as<int>(mattype_scalar_r);

    int dim = std::min(m, n);
    // Fail the case of function-typed matrix
    if(mattype == FUNCTION)
        mattype = -99;
    // Operation for original matrix
    MatProd *op_orig = get_mat_prod_op(A_mat_r, m, n, params_list_r, mattype);
    // Operation for SVD
    MatProd *op;
    if(m > n)
        op = new SVDTallOp(op_orig);
    else
        op = new SVDWideOp(op_orig);

    // Prepare initial residuals
    double *init_resid;
    #include "rands.h"
    if(dim <= rands_len)
    {
        init_resid = rands;
    } else {
        init_resid = new double[dim];
        double *coef_pntr = init_resid;
        for(int i = 0; i < dim / rands_len; i++, coef_pntr += rands_len)
        {
            std::copy(rands, rands + rands_len, coef_pntr);
        }
        std::copy(rands, rands + dim % rands_len, coef_pntr);
    }

    SymEigsSolver<double, LARGEST_ALGE, MatProd> eigs(op, k, ncv);
    eigs.init(init_resid);
    int nconv = eigs.compute(maxitr, tol);
    if(nconv < k)
        Rcpp::warning("only %d singular values converged, less than k = %d", nconv, k);

    nu = std::min(nu, nconv);
    nv = std::min(nv, nconv);

    Eigen::VectorXd evals = eigs.eigenvalues();
    Eigen::MatrixXd evecs = eigs.eigenvectors(std::max(nu, nv));

    Rcpp::NumericVector d(nconv);
    Rcpp::NumericMatrix u(m, nu), v(n, nv);
    int nops = 0;

    // Copy evals to d and take the square root
    std::copy(evals.data(), evals.data() + nconv, d.begin());
    std::transform(d.begin(), d.end(), d.begin(), simple_sqrt);

    // Copy evecs to u or v according to the shape of A
    // If A is tall, copy evecs to v, otherwise copy to u
    if(m > n)
        std::copy(evecs.data(), evecs.data() + nv * n, v.begin());
    else
        std::copy(evecs.data(), evecs.data() + nu * m, u.begin());

    // Calculate the other one
    if(m > n)
    {
        // A = UDV', A'A = VD^2V', AV = UD, ui = A * vi / di
        // evecs has already been copied to v, so we can overwrite evecs
        for(int i = 0; i < nu; i++)
        {
            evecs.col(i) /= d[i];
            op_orig->perform_op(&evecs(0, i), &u(0, i));
            nops++;
        }
    } else {
        // A = UDV', AA' = UD^2U', A'U = VD, vi = A' * ui / di
        // evecs has already been copied to u, so we can overwrite evecs
        for(int i = 0; i < nv; i++)
        {
            evecs.col(i) /= d[i];
            op_orig->perform_tprod(&evecs(0, i), &v(0, i));
            nops++;
        }
    }

    Rcpp::RObject u_ret;
    Rcpp::RObject v_ret;

    if(nu > 0)  u_ret = u;  else  u_ret = R_NilValue;
    if(nv > 0)  v_ret = v;  else  v_ret = R_NilValue;

    if(dim > rands_len)
        delete [] init_resid;

    delete op;
    delete op_orig;

    return Rcpp::List::create(
        Rcpp::Named("d")     = d,
        Rcpp::Named("u")     = u_ret,
        Rcpp::Named("v")     = v_ret,
        Rcpp::Named("niter") = eigs.num_iterations(),
        Rcpp::Named("nops")  = eigs.num_operations() * 2 + nops
    );

    END_RCPP
}
