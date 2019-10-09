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

    MatProd* op = get_mat_prod_op(A_mat_r, n, n, params_list_r, mattype);

    SymEigsSolver<double, LARGEST_MAGN, MatProd> eigs(op, k, ncv);
    eigs.init();
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
                double* ptr = &v(0, i);
                std::transform(ptr, ptr + n, ptr, std::negate<double>());
            }
        }
    }

    Rcpp::RObject u_ret;
    Rcpp::RObject v_ret;

    if(nu > 0)  u_ret = u;  else  u_ret = R_NilValue;
    if(nv > 0)  v_ret = v;  else  v_ret = R_NilValue;

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

    typedef Eigen::Map<const Eigen::VectorXd> MapConstVec;
    typedef Eigen::Map<Eigen::MatrixXd> MapMat;

    Rcpp::List params_svds(params_list_r);

    int m        = as<int>(m_scalar_r);
    int n        = as<int>(n_scalar_r);
    int k        = as<int>(k_scalar_r);
    int nu       = as<int>(nu_scalar_r);
    int nv       = as<int>(nv_scalar_r);
    int ncv      = as<int>(params_svds["ncv"]);
    double tol   = as<double>(params_svds["tol"]);
    int maxitr   = as<int>(params_svds["maxitr"]);
    bool center  = as<bool>(params_svds["center"]);
    bool scale   = as<bool>(params_svds["scale"]);
    int mattype  = as<int>(mattype_scalar_r);

    Rcpp::NumericVector ctr_vec = params_svds["ctr_vec"];
    Rcpp::NumericVector scl_vec = params_svds["scl_vec"];
    MapConstVec ctr_map(ctr_vec.begin(), n);
    MapConstVec scl_map(scl_vec.begin(), n);

    // Operation for original matrix
    MatProd* op_orig = get_mat_prod_op(A_mat_r, m, n, params_list_r, mattype);
    // Operation for SVD
    MatProd* op;
    if(m > n)
        op = new SVDTallOp(op_orig, center, scale, ctr_map, scl_map);
    else
        op = new SVDWideOp(op_orig, center, scale, ctr_map, scl_map);

    SymEigsSolver<double, LARGEST_ALGE, MatProd> eigs(op, k, ncv);
    eigs.init();
    int nconv = eigs.compute(maxitr, tol);
    if(nconv < k)
        Rcpp::warning("only %d singular values converged, less than k = %d", nconv, k);

    nu = std::min(nu, nconv);
    nv = std::min(nv, nconv);

    Eigen::VectorXd evals = eigs.eigenvalues();
    Eigen::MatrixXd evecs = eigs.eigenvectors(std::max(nu, nv));

    Rcpp::NumericVector d(nconv);
    Rcpp::NumericMatrix u(m, nu), v(n, nv);
    MapMat umap = as<MapMat>(u);
    MapMat vmap = as<MapMat>(v);
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
        // B = (A - 1c')S, S = diag(1 / s)
        // B = UDV', B'B = VD^2V', BV = UD, ui = B * vi / di
        // B * x = A * (Sx) - (c'Sx) * 1
        // B * vi / di = A * (vi / s / di) - c'(vi / s / di) * 1
        // evecs has already been copied to v, so we can overwrite evecs
        for(int i = 0; i < nu; i++)
        {
            evecs.col(i).array() /= (d[i] * scl_map).array();
            op_orig->perform_op(&evecs(0, i), &u(0, i));
            umap.col(i).array() -= ctr_map.dot(evecs.col(i));
            nops++;
        }
    } else {
        // B = (A - 1c')S, S = diag(1 / s)
        // B = UDV', BB' = UD^2U', B'U = VD, vi = B' * ui / di
        // B' * x = S(A' * x) - S(c(1'x))
        // B' * ui / di = (A' * (ui / di)) / s - 1'(ui / di) * c / s
        // evecs has already been copied to u, so we can overwrite evecs
        for(int i = 0; i < nv; i++)
        {
            op_orig->perform_tprod(&evecs(0, i), &v(0, i));
            vmap.col(i).array() -= evecs.col(i).sum() * ctr_map.array();
            vmap.col(i).array() /= (d[i] * scl_map.array());
            nops++;
        }
    }

    Rcpp::RObject u_ret;
    Rcpp::RObject v_ret;

    if(nu > 0)  u_ret = u;  else  u_ret = R_NilValue;
    if(nv > 0)  v_ret = v;  else  v_ret = R_NilValue;

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
