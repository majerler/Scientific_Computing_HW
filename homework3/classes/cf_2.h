#ifndef CF_2_H
#define CF_2_H

#include <math.h>

class CF_2
{
public:
    //CF_2();
    virtual double operator()(double x, double y) const = 0;
    // Why use virtual??
    // We are going to overload operator.
    // We are creating a template class.

};

class velocity_X : public CF_2
{
public:
    double operator()(double x, double y) const
    {
        return -1.0*y;
    }
};

class velocity_Y : public CF_2
{
public:
    double operator()(double x, double y) const
    {
        return x;
    }
};

#endif // CF_2_H
