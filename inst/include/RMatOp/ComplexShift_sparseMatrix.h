#ifndef COMPLEXSHIFT_SPARSEMATRIX_H
#define COMPLEXSHIFT_SPARSEMATRIX_H

#include <RcppEigen.h>
#include "ComplexShift.h"

template <int Storage>
class ComplexShift_sparseMatrix: public ComplexShift
{
private:
    typedef Eigen::SparseMatrix<double, Storage> SpMat;
    typedef Eigen::Map<SpMat> MapSpMat;
    typedef Eigen::Map<const Eigen::VectorXd> MapConstVec;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;

    typedef std::complex<double> Complex;
    typedef Eigen::VectorXcd ComplexVector;
    typedef Eigen::SparseMatrix<Complex, Storage> SpCMat;
    typedef Eigen::SparseLU< Eigen::SparseMatrix<Complex, Eigen::ColMajor> > SpLUSolver;

    // Map to Eigen sparse matrix
    MapSpMat      mat;
    const int     n;
    SpLUSolver    solver;
    ComplexVector x_cache;

public:
    ComplexShift_sparseMatrix(SEXP mat_, const int nrow_) :
        mat(Rcpp::as<MapSpMat>(mat_)),
        n(nrow_)
    {}

    int rows() const { return n; }
    int cols() const { return n; }

    void set_shift(double sigmar, double sigmai)
    {
        SpCMat cmat = mat.template cast<Complex>();

        // Create a sparse idendity matrix (1 + 0i on diagonal)
        SpCMat I(n, n);
        I.setIdentity();

        // Sparse LU decomposition
        solver.compute(cmat - Complex(sigmar, sigmai) * I);

        x_cache.resize(n);
        x_cache.setZero();
    }

    // y_out = inv(A - sigma * I) * x_in
    void perform_op(const double* x_in, double* y_out)
    {
        x_cache.real() = MapConstVec(x_in, n);
        MapVec y(y_out, n);
        y.noalias() = solver.solve(x_cache).real();
    }
};

// Operations on "dgCMatrix" class, defined in Matrix package
typedef ComplexShift_sparseMatrix<Eigen::ColMajor> ComplexShift_dgCMatrix;

// Operations on "dgRMatrix" class, defined in Matrix package
typedef ComplexShift_sparseMatrix<Eigen::RowMajor> ComplexShift_dgRMatrix;


#endif // COMPLEXSHIFT_SPARSEMATRIX_H
