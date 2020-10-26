#include "eno_advection.h"
#include <cmath>
//#include <omp.h>


ENO_Advection::ENO_Advection()
{

}


ENO_Advection::ENO_Advection(Grid2D grid_, velocity_X * velx_, velocity_Y * vely_, std::vector <double> sol_tn_)
{
    grid = grid_;
    velx = velx_;
    vely = vely_;
    sol_tn = sol_tn_;
}

std::vector <double> ENO_Advection::get_sol()
{
    return sol_tn;
}

double minmod(double l, double m){
    if (l*m < 0) return 0;
    else if(abs(l) < abs(m)) return l;
    else return m;
    }

void ENO_Advection::one_step_second(double dt)
{
    #pragma omp parallel for
    for (int n = 0; n < grid.get_M() * grid.get_N(); ++n)
    {
        double x = grid.x_from_n(n);
        double y = grid.y_from_n(n);
        int    N = grid.get_N();
        double v_x = (*velx)(x,y);
        double v_y = (*vely)(x,y);
        double dx1, dy1, dx2, dy2;

        if ( v_x < 0)
        {    dx1 = grid.dx_forward(sol_tn, n);
             dx2 = minmod(grid.dx_2(sol_tn, n), grid.dx_2(sol_tn, n + 1));
             }
        if ( v_x >= 0)
        {    dx1 = grid.dx_backward(sol_tn, n);
             dx2 = minmod(grid.dx_2(sol_tn, n), grid.dx_2(sol_tn, n - 1));
             }
        if ( v_y < 0)
        {    dy1 = grid.dy_forward(sol_tn, n);
             dy2 = minmod(grid.dy_2(sol_tn, n), grid.dy_2(sol_tn, n + N));
             }
        if ( v_y >= 0)
        {    dy1 = grid.dy_backward(sol_tn, n);
             dy2 = minmod(grid.dy_2(sol_tn, n), grid.dy_2(sol_tn, n - N ));
             }

        sol_tn[n] = sol_tn[n] - dt * ( v_x * ( (dx1 + dx2 * grid.get_dx()/2) ) )
                              - dt * ( v_y * ( (dy1 + dy2 * grid.get_dy()/2) ) );
    }
}

void ENO_Advection::one_step_central(double dt)
{
    #pragma omp parallel for
    for (int n = 0; n < grid.get_M() * grid.get_N(); ++n)
    {
        double x = grid.x_from_n(n);
        double y = grid.y_from_n(n);
        double v_x = (*velx)(x,y);
        double v_y = (*vely)(x,y);
        double dx1, dy1, dx2, dy2;

        if ( v_x < 0)
        {   dx1 = grid.dx_forward(sol_tn, n);
            dx2 = grid.dx_2(      sol_tn, n); }
        if ( v_x >= 0)
        {   dx1 = grid.dx_backward(sol_tn, n);
            dx2 = grid.dx_2(       sol_tn, n); }
        if ( v_y < 0)
        {   dy1 = grid.dy_forward(sol_tn, n);
            dy2 = grid.dy_2(      sol_tn, n); }
        if ( v_y >= 0)
        {   dy1 = grid.dy_backward(sol_tn, n);
            dy2 = grid.dy_2(       sol_tn, n); }

        sol_tn[n] = sol_tn[n] - dt * v_x * (dx1 + dx2 * grid.get_dx()/2)
                              - dt * v_y * (dy1 + dy2 * grid.get_dy()/2) ;
    }
}

void ENO_Advection::one_step_first(double dt)
{
    #pragma omp parallel for
    for (int n = 0; n < grid.get_M() * grid.get_N(); ++n)
    {
        double x = grid.x_from_n(n);
        double y = grid.y_from_n(n);
        double v_x = (*velx)(x,y);
        double v_y = (*vely)(x,y);
        double dx1, dy1;

        if ( v_x < 0)
        {   dx1 =  grid.dx_forward( sol_tn, n); }
        if ( v_x >= 0)
        {   dx1 =  grid.dx_backward(sol_tn, n); }
        if ( v_y < 0)
        {   dy1 =  grid.dy_forward( sol_tn, n); }
        if ( v_y >= 0)
        {   dy1 =  grid.dy_backward(sol_tn, n); }

        sol_tn[n] = sol_tn[n] - dt * ( v_x * dx1 + v_y * dy1  );
    }
}
