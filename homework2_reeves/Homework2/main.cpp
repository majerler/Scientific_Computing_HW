#include <iostream>
#include <math.h>
#include <cmath>
#include <cf_2.h>
#include <../Lab1/grid2d.h>
#include <eno_advection.h>

using namespace std;

int main()
{
    //grid settings! 1,2,4,8
    int N = 2*64;
    int M = 2*64;
    double xmin = -1.0;
    double xmax = 1.0;
    double ymin = -1.0;
    double ymax = 1.0;

    // Create Grid
    Grid2D grid(N,M,xmin,xmax,ymin,ymax);

    //time settings
    //double dt = grid.get_dx()*grid.get_dx();
    double dt = grid.get_dx() * grid.get_dx();
    double t_f = 4. * acos(0.0);
    cout << "t_f: " << t_f << endl;

    // Create velocity instances
    velocity_X velocity_x;
    velocity_Y velocity_y;

    //Initial solution is a vector form of grid
    //Populate sol vector with initial solution
    vector<double> sol(N*M);

    for (int j = 0; j < N*M; ++j) {
        if ( sqrt( std::pow((grid.x_from_n(j) - 0.5),2) + std::pow(grid.y_from_n(j),2) ) - 0.2 <= 0)
            { sol[j] = 1.0; }
        else
            {sol[j] = 0.0; }
    }

//    vector<double> IC(N*M);

//    for (int j = 0; j < N*M; ++j) {
//        IC[j] = sqrt(std::pow((grid.x_from_n(j) - 0.25),2) + std::pow(grid.y_from_n(j),2))  - 0.2;
//    }



    // Feed grid, velocities, and initital solutions into ENO_Advection class
    ENO_Advection eno_advec(grid, &velocity_x, &velocity_y, sol);
//    ENO_Advection eno_advec(grid, &velocity_x, &velocity_y, IC);

    char name[250];

    double i = 0.0;
    int k = 0;
    while (i < t_f) {

        eno_advec.one_step_central(dt);

        if ( k%100 == 0 ){
            vector<double> sol_comp = eno_advec.get_sol();
            sprintf(name,"/Users/maj/Documents/sol_central/IC_sol%d.vtk", k);
            grid.initialize_VTK_file(name);
            grid.print_VTK_Format(sol_comp,"concentration",name);
            cout << "iterations: " << k << endl;
        }

        i += dt;
        k+= 1;
    }

    cout << "iterations: " << k << endl;
    vector<double> sol_comp = eno_advec.get_sol();

//  Print final solution to vtk file
//    sprintf(name,"/Users/maj/Documents/sols/sol_small_dt_center%d.vtk", N);
//    grid.initialize_VTK_file(name);
//    grid.print_VTK_Format(sol_comp, "value_at_nodes",name);





    //Calculates final solution
    vector<double> sol_final(N*M);

    for (int j = 0; j < N*M; ++j) {
        if ( sqrt( std::pow((grid.x_from_n(j)*cos(t_f) + grid.y_from_n(j)*sin(t_f) - 0.5),2) +
                   std::pow((grid.x_from_n(j)*sin(t_f) + grid.y_from_n(j)*cos(t_f) ),2) ) - 0.2 <= 0)
            { sol_final[j] = 1.0; }
        else
            {sol_final[j] = 0.0; }
    }

    //Error Analysis
    vector<double> dif_between(N*M);
    for (int n = 0; n < N*M; ++n) {
        dif_between[n] = abs(sol_final[n] - sol_comp[n]);
    }

    double error = 0.0;
    for (int n = 0; n < N*M; ++n) {
        error += std::pow(sol_final[n] - sol_comp[n], 2);
    }


//        sprintf(name,"/Users/maj/Documents/sols/diff%d.vtk", N);
//        grid.initialize_VTK_file(name);
//        grid.print_VTK_Format(dif_between, "value_at_nodes",name);

    double norm_error = sqrt(error) / double (N*M);
    cout << "error: " << norm_error << endl;


    return 0;

}

