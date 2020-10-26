#ifndef GODUNOV_H
#define GODUNOV_H

#include <../classes/grid2d.h>
#include <../classes/cf_2.h>
#include <../classes/math_tools.h>
#include <vector>

class godunov
{
private:
    Grid2D grid;
    std::vector <double> sol;
    std::vector <double> og;
    velocity_X *velx;
    velocity_Y *vely;
public:
    godunov();
    godunov(Grid2D grid_, velocity_X * velx_, velocity_Y * vely_,  std::vector <double> solution_, std::vector <double> original_);
    void upwind(double dt);
    std::vector <double> get_sol();
};

#endif // GODUNOV_H
