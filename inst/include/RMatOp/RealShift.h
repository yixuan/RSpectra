#ifndef REALSHIFT_H
#define REALSHIFT_H

class RealShift
{
public:
    virtual int rows() = 0;
    virtual int cols() = 0;

    virtual void set_shift(double sigma) = 0;

    // y_out = inv(A - sigma * I) * x_in
    virtual void perform_op(double *x_in, double *y_out) = 0;

    virtual ~RealShift() {}
};


#endif // REALSHIFT_H
