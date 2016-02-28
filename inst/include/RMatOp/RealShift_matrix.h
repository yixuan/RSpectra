#ifndef REALSHIFT_MATRIX_H
#define REALSHIFT_MATRIX_H

#include <RcppEigen.h>
#include "RealShift.h"

class RealShift_matrix: public RealShift
{
private:
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::Map<Eigen::MatrixXd> MapMat;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;
    typedef Eigen::PartialPivLU<Eigen::MatrixXd> LUSolver;

    MapMat mat;
    const int n;
    LUSolver solver;

public:
    RealShift_matrix(SEXP mat_, const int nrow_) :
        mat(REAL(mat_), nrow_, nrow_),
        n(nrow_)
    {}

    int rows() { return n; }
    int cols() { return n; }

    void set_shift(double sigma)
    {
        solver.compute(mat - sigma * Matrix::Identity(n, n));
    }

    // y_out = inv(A - sigma * I) * x_in
    void perform_op(double *x_in, double *y_out)
    {
        MapVec x(x_in, n);
        MapVec y(y_out, n);
        y.noalias() = solver.solve(x);
    }
};


#endif // REALSHIFT_MATRIX_H
