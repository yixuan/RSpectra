#ifndef REALSHIFT_SYM_SPARSEMATRIX_H
#define REALSHIFT_SYM_SPARSEMATRIX_H

#include <RcppEigen.h>
#include <Rdefines.h>  // for R macros
#include "RealShift.h"

template <int Storage>
class RealShift_sym_sparseMatrix: public RealShift
{
private:
    typedef Eigen::SparseMatrix<double, Storage> SpMat;
    typedef Eigen::MappedSparseMatrix<double, Storage> MapSpMat;
    typedef Eigen::Map<const Eigen::VectorXd> MapConstVec;
    typedef Eigen::Map<Eigen::VectorXd> MapVec;
    typedef Eigen::SimplicialLDLT< Eigen::SparseMatrix<double, Eigen::ColMajor> > SpLDLSolver;

    // Map to Eigen sparse matrix
    MapSpMat    mat;
    const int   n;
    const char  uplo;
    SpLDLSolver solver;

public:
    RealShift_sym_sparseMatrix(SEXP mat_, const int nrow_) :
        mat(Rcpp::as<MapSpMat>(mat_)),
        n(nrow_),
        uplo(Rcpp::as<std::string>(GET_SLOT(mat_, Rf_install("uplo"))))
    {}

    int rows() const { return n; }
    int cols() const { return n; }

    void set_shift(double sigma)
    {
        solver.setShift(-sigma);

        if(uplo == 'L')
            solver.compute(mat.template selfadjointView<Eigen::Lower>());
        else
            solver.compute(mat.template selfadjointView<Eigen::Upper>());
    }

    // y_out = inv(A - sigma * I) * x_in
    void perform_op(const double* x_in, double* y_out)
    {
        MapConstVec x(x_in, n);
        MapVec y(y_out, n);
        y.noalias() = solver.solve(x);
    }
};

// Operations on "dgCMatrix" class, defined in Matrix package
typedef RealShift_sym_sparseMatrix<Eigen::ColMajor> RealShift_dsCMatrix;

// Operations on "dgRMatrix" class, defined in Matrix package
typedef RealShift_sym_sparseMatrix<Eigen::RowMajor> RealShift_dsRMatrix;


#endif // REALSHIFT_SYM_SPARSEMATRIX_H
