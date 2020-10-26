#ifndef SL_METHOD_H
#define SL_METHOD_H

#include <../classes/grid2d.h>
#include <../classes/cf_2.h>
#include <../classes/math_tools.h>
//#include <../../../hw3/omp.h>
#include <omp.h>
#include <vector>

class SL_method
{
private:
    Grid2D grid;
    std::vector <double> sol;
    std::vector <double> og;
    velocity_X *velx ;
    velocity_Y *vely ;
public:
    SL_method();
    SL_method(Grid2D grid_, velocity_X * velx_, velocity_Y * vely_,  std::vector <double> solution_, std::vector <double> og_);
    double clipclip(double in, int j);
    std::vector<double> trajectory_interpolation(int n, double dt);
    void one_step(double dt);
    void godunov(double dt);
    std::vector <double> get_sol();
    void reset_og(std::vector <double> og_);
};

#endif // SL_METHOD_H
