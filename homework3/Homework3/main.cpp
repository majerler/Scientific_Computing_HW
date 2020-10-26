
#include <iostream>
#include <../classes/grid2d.h>
#include <../classes/math_tools.h>
#include <../classes/sl_method.h>
#include <../classes/godunov.h>

#include <vector>
#include <math.h>
#include <cmath>

using namespace std;

int main()
{
    cout << "Hello World!" << endl;

    // Set up for grid
    double xmin = -1;
    double xmax = 1;
    double ymin = -1;
    double ymax = 1;

    int N = 1*64;
    int M = 1*64;

    // create grid and velocity objects
    Grid2D grid(N, M, xmin, xmax, ymin, ymax);
    velocity_X velocity_x;
    velocity_Y velocity_y;

    // create initial solution
    vector <double> initial_sol(N*M);
    vector <double> pert_sol(N*M);
    for (int n = 0; n < N*M; ++n)
    {
        initial_sol[n] = sqrt( pow(grid.x_from_n(n) - 0.25 ,2) + pow(grid.y_from_n(n),2) ) - 0.2 ;
        pert_sol[n] = 2.0 * initial_sol[n]  ;
    }


//    vector <double> initial_sol(N*M);
//    for (int n = 0; n < N*M; ++n)
//    {
//        double y = grid.y_from_n(n);
//        double x = grid.x_from_n(n);
//        if ( y > -0.25 & y < 2*x + 0.25 & y < -2*x + 0.25 ){
//            initial_sol[n] = 1;
//        }
//        else{
//            initial_sol[n] = 0;
//        }

//    }

    // Input everything into semi_lagrangian and gudonov method
    SL_method semi_lagrange(grid, &velocity_x, &velocity_y, initial_sol, initial_sol);
//    SL_method semi_lagrange(grid, &velocity_x, &velocity_y, pert_sol, initial_sol);

    // set up for stepping through time
    double t_fin = 4 * acos(0);
    double dt =  grid.get_dx() / 15.;
    double dt2 =  grid.get_dx() / 2.0;
    int num_iter = t_fin / dt;
    double t_f = dt * (double) num_iter;
    cout << "final time experiment " << t_f << " true final time " << t_fin <<endl;

    //  ##################################### TEST SEMI - LAGRANGIAN #####################################
/*
        char name[250];
        for (int i = 0; i < num_iter; ++i) {

            semi_lagrange.one_step(dt);

            if ( i % 100 ){
                vector<double> final_sol = semi_lagrange.get_sol();
                sprintf(name,"/Users/maj/Desktop/hw3/ec/movie/iter_sol%d.vtk", i);
                grid.initialize_VTK_file(name);
                grid.print_VTK_Format(final_sol,"concentration",name);
            }
        }

        vector <double> final_sol = semi_lagrange.get_sol();

        char name_f[250];
        sprintf(name_f,"/Users/maj/Desktop/hw3/ec/final%d.vtk", N);
        grid.initialize_VTK_file(name_f);
        grid.print_VTK_Format(final_sol, "value_at_nodes",name_f);

        char name_i[250];
        sprintf(name_i,"/Users/maj/Desktop/hw3/ec/initial%d.vtk", N);
        grid.initialize_VTK_file(name_i);
        grid.print_VTK_Format(initial_sol, "value_at_nodes",name_i);
*/


    //  ##################################### TEST SEMI - LAGRANGIAN #####################################

        char name[250];
        for (int i = 0; i < num_iter; ++i) {

            semi_lagrange.one_step(dt);

            if ( i%100 == 0 ){
                vector<double> final_sol = semi_lagrange.get_sol();
                sprintf(name,"/Users/maj/Desktop/hw3/sl/movie/iter_sol%d.vtk", i);
                grid.initialize_VTK_file(name);
                grid.print_VTK_Format(final_sol,"concentration",name);
            }
        }

        vector <double> final_sol = semi_lagrange.get_sol();

        //Calculates final solution
        vector<double> true_sol(N*M);
        for (int j = 0; j < N*M; ++j) {
            true_sol[j] = sqrt( pow(grid.x_from_n(j)*cos(t_fin) + grid.y_from_n(j)*sin(t_fin) - 0.25,2) +
                                 pow(grid.y_from_n(j)*cos(t_fin) - grid.x_from_n(j)*sin(t_fin) ,2) ) - 0.2;
        }

        vector<double> dif_between(N*M);
        for (int n = 0; n < N*M; ++n) {
            dif_between[n] = abs(true_sol[n] - final_sol[n]); }


        double error = 0.0;
        for (int n = 0; n < N*M; ++n) {
            error += pow(true_sol[n] - final_sol[n], 2); }

        double norm_error = sqrt(error) / (double) (N*M);
        cout << "error: " << norm_error << endl;

        char name_i[250];
        sprintf(name_i,"/Users/maj/Desktop/hw3/sl/initial%d.vtk", N);
        grid.initialize_VTK_file(name_i);
        grid.print_VTK_Format(initial_sol, "value_at_nodes",name_i);

        char name_t[250];
        sprintf(name_t,"/Users/maj/Desktop/hw3/sl/true%d.vtk", N);
        grid.initialize_VTK_file(name_t);
        grid.print_VTK_Format(true_sol, "value_at_nodes",name_t);

        char name_f[250];
        sprintf(name_f,"/Users/maj/Desktop/hw3/sl/final%d.vtk", N);
        grid.initialize_VTK_file(name_f);
        grid.print_VTK_Format(final_sol, "value_at_nodes",name_f);

        char name_d[250];
        sprintf(name_d,"/Users/maj/Desktop/hw3/sl/diff%d.vtk", N);
        grid.initialize_VTK_file(name_d);
        grid.print_VTK_Format(dif_between, "value_at_nodes",name_d);

        double error2 = 0;
        for (int i = 0; i < N*M; ++i) {
            if (dif_between[i] >= error2){
                error2 = dif_between[i];
            }
        }
        cout << "infinity error: " << error2 << endl;



    //  ##################################### TEST GODUNOV #####################################
    /*
        char name[250];
        for (int i = 0; i < 500; ++i) {
            semi_lagrange.godunov(dt2) ;

            vector<double> final_sol = semi_lagrange.get_sol();
            sprintf(name,"/Users/maj/Desktop/hw3/reinitial/movie/iter_sol%d.vtk", i);
            grid.initialize_VTK_file(name);
            grid.print_VTK_Format(final_sol,"concentration",name);

        }


        vector <double> final_sol = semi_lagrange.get_sol();

        vector<double> dif_between(N*M);
        for (int n = 0; n < N*M; ++n) {
            dif_between[n] = abs(initial_sol[n] - final_sol[n]); }


        char name_i[250];
        sprintf(name_i,"/Users/maj/Desktop/hw3/reinitial/initial%d.vtk", N);
        grid.initialize_VTK_file(name_i);
        grid.print_VTK_Format(initial_sol, "value_at_nodes",name_i);

        char name_f[250];
        sprintf(name_f,"/Users/maj/Desktop/hw3/reinitial/final%d.vtk", N);
        grid.initialize_VTK_file(name_f);
        grid.print_VTK_Format(final_sol, "value_at_nodes",name_f);

        char name_d[250];
        sprintf(name_d,"/Users/maj/Desktop/hw3/reinitial/diff%d.vtk", N);
        grid.initialize_VTK_file(name_d);
        grid.print_VTK_Format(dif_between, "value_at_nodes",name_d);

    */

//  ##################################### TEST SL w/ GODUNOV #####################################
/*
    char name[250];
    for (int i = 0; i < num_iter; ++i) {

        semi_lagrange.one_step(dt);

        for (int j = 0; j < 5; ++j) {
            semi_lagrange.godunov(dt2);
        }

        if ( N == 64 & dt == grid.get_dx() / 10. ){
            vector<double> final_sol = semi_lagrange.get_sol();
            sprintf(name,"/Users/maj/Desktop/hw3/joint/movie/iter_sol%d.vtk", i);
            grid.initialize_VTK_file(name);
            grid.print_VTK_Format(final_sol,"concentration",name);
        }

//        phi_0 = semi_lagrange.get_sol();
//        semi_lagrange.reset_og(phi_0);
    }

    vector <double> final_sol = semi_lagrange.get_sol();

    //Calculates final solution
    vector<double> true_sol(N*M);
    for (int j = 0; j < N*M; ++j) {
        true_sol[j] = sqrt( pow(grid.x_from_n(j)*cos(t_fin) + grid.y_from_n(j)*sin(t_fin) - 0.25,2) +
                            pow(grid.y_from_n(j)*cos(t_fin) - grid.x_from_n(j)*sin(t_fin) ,2) ) - 0.2;
    }

    vector<double> dif_between(N*M);
    for (int n = 0; n < N*M; ++n) {
        dif_between[n] = abs(true_sol[n] - final_sol[n]); }


    double error = 0.0;
    for (int n = 0; n < N*M; ++n) {
        error += pow(true_sol[n] - final_sol[n], 2); }

    double norm_error = sqrt(error) / (double) (N*M);
    cout << "error: " << norm_error << endl;

    char name_i[250];
    sprintf(name_i,"/Users/maj/Desktop/hw3/joint/initial%d.vtk", N);
    grid.initialize_VTK_file(name_i);
    grid.print_VTK_Format(initial_sol, "value_at_nodes",name_i);

    char name_t[250];
    sprintf(name_t,"/Users/maj/Desktop/hw3/joint/true%d.vtk", N);
    grid.initialize_VTK_file(name_t);
    grid.print_VTK_Format(true_sol, "value_at_nodes",name_t);

    char name_f[250];
    sprintf(name_f,"/Users/maj/Desktop/hw3/joint/final%d.vtk", N);
    grid.initialize_VTK_file(name_f);
    grid.print_VTK_Format(final_sol, "value_at_nodes",name_f);

    char name_d[250];
    sprintf(name_d,"/Users/maj/Desktop/hw3/joint/diff%d.vtk", N);
    grid.initialize_VTK_file(name_d);
    grid.print_VTK_Format(dif_between, "value_at_nodes",name_d);
*/



    return 0;
}
