#include "sl_method.h"
#include <math.h>
//#include <../classes/grid2d.h>
//#include <../classes/cf_2.h>

SL_method::SL_method()
{

}

SL_method::SL_method(Grid2D grid_,  velocity_X * velx_, velocity_Y * vely_, std::vector <double> solution_, std::vector <double> og_)
{
    grid = grid_;
    sol = solution_;
    velx = velx_;
    vely = vely_;
    og = og_;

}

void SL_method::reset_og(std::vector <double> og_)
{
    og = og_;
}

double SL_method::clipclip(double in, int j){

    double clip = in;

    if (j == 0){
        double xmin = grid.get_xmin();
        double xmax = grid.get_xmax();

        clip = std::max(xmin, clip);
        clip = std::min(xmax, clip);
    }

    if (j == 1){
        double ymin = grid.get_ymin();
        double ymax = grid.get_ymax();

        clip = std::max(ymin, clip);
        clip = std::min(ymax, clip);
    }

    return clip;
}

std::vector<double> SL_method::trajectory_interpolation(int n, double dt)
{
    std::vector<double> xy_d(2);

    double x = grid.x_from_n(n);
    double y = grid.y_from_n(n);

    double x_star = x - dt * (*velx)(x,y) / 2;
    double y_star = y - dt * (*vely)(x,y) / 2;

//    x_star = clipclip(x_star, 0);
//    y_star = clipclip(y_star, 1);

    xy_d[0] =       x - dt * (*velx)(x_star,y_star);
    xy_d[1] =       y - dt * (*vely)(x_star,y_star);

//    xy_d[0] = clipclip(xy_d[0], 0);
//    xy_d[1] = clipclip(xy_d[1], 1);

    return xy_d;

}

void SL_method::one_step(double dt)
{
    #pragma omp parallel for
    for (int n = 0; n < grid.get_N() * grid.get_M(); ++n) {
        std::vector<double> xy_depart = trajectory_interpolation(n, dt);
        sol[n] = eno_interpolation(grid, sol, xy_depart[0], xy_depart[1]);
        og[n] = sol[n];
    }

}

void SL_method::godunov(double dt)
{
    #pragma omp parallel for
    for (int n = 0; n < grid.get_M() * grid.get_N(); ++n)
    {
        double dx1, dy1;

        double dxf = grid.dx_forward( sol, n);
        double dxb = grid.dx_backward( sol, n);
        double dyf = grid.dy_forward( sol, n);
        double dyb = grid.dy_backward( sol, n);
        double s;


//        if ( og[n]     >  1 * grid.get_dx() ){
//            s = 1.;         }
//        else if (og[n] < -1 * grid.get_dx() ) {
//            s = -1.;        }
//        else{
//            s = 0.;         }

        if ( og[n]     >  0.25 * grid.get_dx() ){
            s = 1.;         }
        else if (og[n] < - 0.25 * grid.get_dx() ) {
            s = -1.;        }
        else{
            s = 0.;         }


        if ( s*dxb <= 0 & s*dxf <= 0 ){
            dx1 = dxf ;             }

        if ( s*dxb >= 0 & s*dxf >= 0 ){
            dx1 = dxb ;             }

        if ( s*dxb <= 0 & s*dxf >= 0 ){
            dx1 = 0.;               }

        if ( s*dxb >= 0 & s*dxf <= 0 ){
            if( abs(dxb) >= abs(dxf)){
                dx1 = dxb ;         }
            else if ( abs(dxb) <= abs(dxf)){
                dx1 = dxf ;         }
        }

        if ( s*dyb <= 0 & s*dyf <= 0 ){
            dy1 = dyf ;             }

        if ( s*dyb >= 0 & s*dyf >= 0 ){
            dy1 = dyb ;             }

        if ( s*dyb <= 0 & s*dyf >= 0 ){
            dy1 = 0. ;              }

        if ( s*dyb >= 0 & s*dyf <= 0 ){
            if( abs(dyb) >= abs(dyf)){
                dy1 = dyb ;         }
            else if ( abs(dyb) <= abs(dyf)){
                dy1 = dyf ;         }
        }

        sol[n] = sol[n] + dt*s* ( 1 - sqrt(dx1*dx1 + dy1*dy1) );

    }
}


std::vector <double> SL_method::get_sol()
{
    return sol;
}


