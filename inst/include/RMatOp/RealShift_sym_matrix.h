#ifndef REALSHIFT_SYM_MATRIX_H
#define REALSHIFT_SYM_MATRIX_H

#include <RcppEigen.h>

class RealShift_sym_matrix: public RealShift
{
private:
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::VectorXd Vector;
    typedef Eigen::Map<Eigen::MatrixXd> MapMat;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;
    typedef Eigen::LDLT<Eigen::MatrixXd> LDLTSolver;

    MapMat mat;
    const int n;
    const char uplo;
    LDLTSolver solver;

public:
    RealShift_sym_matrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        mat(REAL(mat_), nrow_, nrow_),
        n(nrow_),
        uplo(uplo_)
    {}

    int rows() { return n; }
    int cols() { return n; }

    void set_shift(double sigma)
    {
        // Backup the diagonal elements
        Vector diag = mat.diagonal();
        mat.diagonal().array() -= sigma;

        if(uplo == 'L')
            solver.compute(mat.selfadjointView<Eigen::Lower>());
        else
            solver.compute(mat.selfadjointView<Eigen::Upper>());

        mat.diagonal() = diag;
    }

    // y_out = inv(A - sigma * I) * x_in
    void perform_op(double *x_in, double *y_out)
    {
        MapVec x(x_in, n);
        MapVec y(y_out, n);
        y.noalias() = solver.solve(x);
    }
};


#endif // REALSHIFT_SYM_MATRIX_H
