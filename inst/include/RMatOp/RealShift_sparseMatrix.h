#ifndef REALSHIFT_SPARSEMATRIX_H
#define REALSHIFT_SPARSEMATRIX_H

#include <RcppEigen.h>

template <int Storage>
class RealShift_sparseMatrix: public RealShift
{
private:
    typedef Eigen::SparseMatrix<double, Storage> SpMat;
    typedef Eigen::MappedSparseMatrix<double, Storage> MapSpMat;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;
    typedef Eigen::SparseLU< Eigen::SparseMatrix<double, Eigen::ColMajor> > SpLUSolver;

    // Map to Eigen sparse matrix
    MapSpMat mat;
    const int n;
    SpLUSolver solver;

public:
    RealShift_sparseMatrix(SEXP mat_, const int nrow_) :
        mat(Rcpp::as<MapSpMat>(mat_)),
        n(nrow_)
    {}

    int rows() { return n; }
    int cols() { return n; }

    void set_shift(double sigma)
    {
        // Create a sparse idendity matrix
        SpMat I(n, n);
        I.setIdentity();

        solver.compute(mat - sigma * I);
    }

    // y_out = inv(A - sigma * I) * x_in
    void perform_op(double *x_in, double *y_out)
    {
        MapVec x(x_in, n);
        MapVec y(y_out, n);
        y.noalias() = solver.solve(x);
    }
};

// Operations on "dgCMatrix" class, defined in Matrix package
typedef RealShift_sparseMatrix<Eigen::ColMajor> RealShift_dgCMatrix;

// Operations on "dgRMatrix" class, defined in Matrix package
typedef RealShift_sparseMatrix<Eigen::RowMajor> RealShift_dgRMatrix;


#endif // REALSHIFT_SPARSEMATRIX_H
