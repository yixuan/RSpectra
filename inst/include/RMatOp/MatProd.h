#ifndef MATPROD_H
#define MATPROD_H

class MatProd
{
public:
    virtual int rows() const = 0;
    virtual int cols() const = 0;

    // y_out = A * x_in
    virtual void perform_op(const double* x_in, double* y_out) = 0;

    // y_out = A' * x_in
    virtual void perform_tprod(const double* x_in, double* y_out) = 0;

    virtual ~MatProd() {}
};

#endif // MATPROD_H
