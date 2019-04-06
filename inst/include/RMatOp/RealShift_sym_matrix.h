#ifndef REALSHIFT_SYM_MATRIX_H
#define REALSHIFT_SYM_MATRIX_H

#include <RcppEigen.h>
#include <R_ext/Lapack.h>
#include "RealShift.h"

class RealShift_sym_matrix: public RealShift
{
private:
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::VectorXd Vector;
    typedef Eigen::VectorXi IntVector;
    typedef Eigen::Map<Eigen::MatrixXd> MapMat;
    typedef Eigen::Map<const Matrix> MapConstMat;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;

    const int  n;
    const char uplo;
    Matrix     fac;   // store the LDLT factorization result
    IntVector  perm;  // store the LDLT factorization result

public:
    RealShift_sym_matrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        n(nrow_), uplo(uplo_), fac(nrow_, nrow_), perm(nrow_)
    {
        fac.noalias() = MapConstMat(REAL(mat_), nrow_, nrow_);
    }

    int rows() const { return n; }
    int cols() const { return n; }

    // Use LAPACK to do the LDLT factorization
    void set_shift(double sigma)
    {
        fac.diagonal().array() -= sigma;

        int lwork = -1, info;
        double blocksize;

        // Inquire the optimal block size
        F77_CALL(dsytrf)(&uplo, &n,
                 fac.data(), &n,
                 perm.data(),
                 &blocksize, &lwork, &info);

        if(info != 0)
            Rcpp::stop("RealShift_sym_matrix: factorization failed with the given shift");

        lwork = int(blocksize);
        std::vector<double> work;
        work.resize(lwork);

        // LDLT factorization
        F77_CALL(dsytrf)(&uplo, &n,
                 fac.data(), &n,
                 perm.data(),
                 &work[0], &lwork, &info);

        if(info != 0)
            Rcpp::stop("RealShift_sym_matrix: factorization failed with the given shift");
    }

    // y_out = inv(A - sigma * I) * x_in
    void perform_op(const double* x_in, double* y_out)
    {
        std::copy(x_in, x_in + n, y_out);

        const int nrhs = 1;
        int info;
        F77_CALL(dsytrs)(&uplo, &n, &nrhs,
                 fac.data(), &n,
                 perm.data(), y_out, &n, &info);

        if(info != 0)
            Rcpp::stop("RealShift_sym_matrix: input vector has illegal values");
    }
};


#endif // REALSHIFT_SYM_MATRIX_H
