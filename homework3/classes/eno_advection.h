#ifndef ENO_ADVECTION_H
#define ENO_ADVECTION_H

#include <../classes/grid2d.h>
#include <../classes/cf_2.h>
#include <vector>


class ENO_Advection
{
private:
    Grid2D grid;
    std::vector <double> sol_tn;
    // These will point to something eventually
//    velocity_X *velx ;
//    velocity_Y *vely ;
    velocity_X *velx ;
    velocity_Y *vely ;
public:
    ENO_Advection();
    std::vector <double> get_sol();
    //ENO_Advection(Grid2D grid_, velocity_X * velx_, velocity_Y * vely_, std::vector <double> sol_tn_);
    ENO_Advection(Grid2D grid_,velocity_X *velx_,velocity_Y *vely_, std::vector <double> sol_tn_);
    void one_step_second(double dt);
    void one_step_central(double dt);
    void one_step_first(double dt);
};



#endif // ENO_ADVECTION_H
