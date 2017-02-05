#ifndef MATOPS_C_H
#define MATOPS_C_H

#include <SpectraC.h>

class CMatProd: public MatProd
{
private:
    mat_op    op;
    const int n;
    void*     data;
public:
    CMatProd(mat_op op_, int n_, void *data_) :
        op(op_),
        n(n_),
        data(data_)
    {}
    int rows() const { return n; }
    int cols() const { return n; }
    void perform_op(const double* x_in, double* y_out) { op(x_in, y_out, n, data); }
    void perform_tprod(const double* x_in, double* y_out) {}
};

class CRealShift: public RealShift
{
private:
    mat_op    op;
    const int n;
    void*     data;
public:
    CRealShift(mat_op op_, int n_, void *data_) :
        op(op_),
        n(n_),
        data(data_)
    {}
    int rows() const { return n; }
    int cols() const { return n; }
    void set_shift(double sigma) {}
    void perform_op(const double* x_in, double* y_out) { op(x_in, y_out, n, data); }
};

class CComplexShift: public ComplexShift
{
private:
    mat_op    op;
    const int n;
    void*     data;
public:
    CComplexShift(mat_op op_, int n_, void *data_) :
        op(op_),
        n(n_),
        data(data_)
    {}
    int rows() const { return n; }
    int cols() const { return n; }
    void set_shift(double sigmar, double sigmai) {}
    void perform_op(const double* x_in, double* y_out) { op(x_in, y_out, n, data); }
};


#endif // MATOPS_C_H
