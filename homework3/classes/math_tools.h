#include <../classes/grid2d.h>
#include <vector>



double bilinear_interpolation(Grid2D & grid, std::vector <double> & func, double x, double y);
double eno_interpolation(Grid2D & grid, std::vector <double> & func, double x, double y);
double minmod(double l, double m);
