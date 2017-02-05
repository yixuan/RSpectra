#ifndef COMPLEXSHIFT_MATRIX_H
#define COMPLEXSHIFT_MATRIX_H

#include <RcppEigen.h>
#include "ComplexShift.h"

class ComplexShift_matrix: public ComplexShift
{
private:
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::Map<const Eigen::MatrixXd> MapConstMat;
    typedef Eigen::Map<const Eigen::VectorXd> MapConstVec;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;

    typedef std::complex<double> Complex;
    typedef Eigen::MatrixXcd ComplexMatrix;
    typedef Eigen::VectorXcd ComplexVector;
    typedef Eigen::PartialPivLU<ComplexMatrix> ComplexSolver;

    MapConstMat   mat;
    const int     n;
    ComplexSolver solver;
    ComplexVector x_cache;

public:
    ComplexShift_matrix(SEXP mat_, const int nrow_) :
        mat(REAL(mat_), nrow_, nrow_),
        n(nrow_)
    {}

    int rows() const { return n; }
    int cols() const { return n; }

    void set_shift(double sigmar, double sigmai)
    {
        ComplexMatrix cmat = mat.cast<Complex>();
        cmat.diagonal().array() -= Complex(sigmar, sigmai);
        solver.compute(cmat);
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


#endif // COMPLEXSHIFT_MATRIX_H
