#ifndef NET_TOOLS_H
#define NET_TOOLS_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <net_tools.h>

#endif // NET_TOOLS_H

double relu(double x);
double sigmoid(double x);
double relu_dx(double x);
double sigmoid_dx(double x);
double cost(std::vector<double> out, std::vector<double> in);
