#include "net_tools.h"
#include <math.h>

double relu(double x)
{
    return fmax(0., x);
}

double sigmoid(double x)
{
    return 1 / ( 1 + exp(-x) );
}

double relu_dx(double x)
{
    if (x > 0)
    {
        return 1.0 ;
    }
    else
    {
    return 0.0 ;
    }
}

double sigmoid_dx(double x)
{
    return sigmoid(x) * ( 1 - sigmoid(x) );
}

double cost(std::vector<double> out, std::vector<double> in)
{
    double error = 0.0;

    for (int j = 0; j < out.size(); ++j) {
        error += 0.5 * ( out[j] - in[j] ) * ( out[j] - in[j] );
    }

    return error;
}
