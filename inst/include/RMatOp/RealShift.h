#ifndef REALSHIFT_H
#define REALSHIFT_H

class RealShift
{
public:
    virtual int rows() const = 0;
    virtual int cols() const = 0;

    virtual void set_shift(double sigma) = 0;

    // y_out = inv(A - sigma * I) * x_in
    virtual void perform_op(const double* x_in, double* y_out) = 0;

    virtual ~RealShift() {}
};


#endif // REALSHIFT_H
