#ifndef SVDOP_H
#define SVDOP_H

#include <RcppEigen.h>
#include <algorithm>
#include "MatProd.h"

class SVDTallOp: public MatProd
{
private:
    typedef Eigen::VectorXd Vector;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;

    MatProd*    op;
    const int   nrow;
    const int   ncol;
    const int   dim;

    const bool  center;
    const bool  scale;
    MapConstVec ctr_vec;
    MapConstVec scl_vec;

    Vector      workm;
    Vector      workn;
public:
    SVDTallOp(MatProd* op_, bool center_, bool scale_,
              const MapConstVec& ctr_vec_, const MapConstVec& scl_vec_) :
        op(op_), nrow(op->rows()), ncol(op->cols()), dim(ncol),
        center(center_), scale(scale_),
        ctr_vec(ctr_vec_.data(), ctr_vec_.size()),
        scl_vec(scl_vec_.data(), scl_vec_.size()),
        workm(nrow), workn(ncol)
    {}

    // Number of rows and columns of the operator B'B, not of B itself
    int rows() const { return dim; }
    int cols() const { return dim; }

    // y_out = B'B * x_in
    // B = (A - 1c')S, c = ctr_vec, S = diag(1 / scl_vec)
    // Bv = A * (Sv) - (c'Sv) * 1
    // B'v = S(A'v) - S((1'v)c)
    void perform_op(const double* x_in, double* y_out)
    {
        // No centering or scaling
        if(!(center || scale))
        {
            op->perform_op   (x_in, workm.data());
            op->perform_tprod(workm.data(), y_out);
        } else {
            // Sv = v / s
            workn.noalias() = MapConstVec(x_in, ncol).cwiseQuotient(scl_vec);
            // A * (Sv)
            op->perform_op(workn.data(), workm.data());
            // A * (Sv) - (c'Sv) * 1
            workm.array() -= ctr_vec.dot(workn);

            // A'v
            op->perform_tprod(workm.data(), workn.data());
            // A'v - (1'v)c
            workn.noalias() -= workm.sum() * ctr_vec;
            // S(A'v) - S((1'v)c)
            MapVec(y_out, ncol).array() = workn.array() / scl_vec.array();
        }
    }

    void perform_tprod(const double* x_in, double* y_out)
    {
        perform_op(x_in, y_out);
    }
};

class SVDWideOp: public MatProd
{
private:
    typedef Eigen::VectorXd Vector;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;

    MatProd*    op;
    const int   nrow;
    const int   ncol;
    const int   dim;

    const bool  center;
    const bool  scale;
    MapConstVec ctr_vec;
    MapConstVec scl_vec;

    Vector      workm;
    Vector      workn;
public:
    SVDWideOp(MatProd* op_, bool center_, bool scale_,
              const MapConstVec& ctr_vec_, const MapConstVec& scl_vec_) :
        op(op_), nrow(op->rows()), ncol(op->cols()), dim(nrow),
        center(center_), scale(scale_),
        ctr_vec(ctr_vec_.data(), ctr_vec_.size()),
        scl_vec(scl_vec_.data(), scl_vec_.size()),
        workm(nrow), workn(ncol)
    {}

    // Number of rows and columns of the operator BB', not of B itself
    int rows() const { return dim; }
    int cols() const { return dim; }

    // y_out = BB' * x_in
    // B = (A - 1c')S, c = ctr_vec, S = diag(1 / scl_vec)
    // B'v = S(A'v) - S((1'v)c)
    // Bv = A * (Sv) - (c'Sv) * 1
    void perform_op(const double* x_in, double* y_out)
    {
        if(!(center || scale))
        {
            op->perform_tprod(x_in, workn.data());
            op->perform_op(workn.data(), y_out);
        } else {
            // A'v
            op->perform_tprod(x_in, workn.data());
            // A'v - (1'v)c
            workn.noalias() -= MapConstVec(x_in, nrow).sum() * ctr_vec;
            // S( S(A'v) - S((1'v)c) ) = S^2(A'v - (1'v)c)
            workn.array() /= scl_vec.array().square();

            // A * (Sv)
            op->perform_op(workn.data(), y_out);
            // A * (Sv) - (c'Sv) * 1
            MapVec(y_out, nrow).array() -= ctr_vec.dot(workn);
        }
    }

    void perform_tprod(const double* x_in, double* y_out)
    {
        perform_op(x_in, y_out);
    }
};

#endif // SVDOP_H
