#ifndef MATPROD_H
#define MATPROD_H

class MatProd
{
public:
    virtual int rows() = 0;
    virtual int cols() = 0;

    // y_out = A * x_in
    virtual void perform_op(double *x_in, double *y_out) = 0;

    // y_out = A' * x_in
    virtual void perform_tprod(double *x_in, double *y_out) = 0;

    virtual ~MatProd() {}
};

#endif // MATPROD_H
