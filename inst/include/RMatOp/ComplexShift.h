#ifndef COMPLEXSHIFT_H
#define COMPLEXSHIFT_H

class ComplexShift
{
public:
    typedef double Scalar;

    virtual int rows() const = 0;
    virtual int cols() const = 0;

    virtual void set_shift(double sigmar, double sigmai) = 0;

    // y_out = Re( inv(A - sigma * I) * x_in )
    virtual void perform_op(const double* x_in, double* y_out) const = 0;

    virtual ~ComplexShift() {}
};


#endif // COMPLEXSHIFT_H
