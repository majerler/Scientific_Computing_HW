#ifndef GRID2D_H
#define GRID2D_H

#include <iostream>
#include <vector>


class Grid2D
{
private:
    int M;
    double xmin, xmax, ymin, ymax;
    double dx, dy;
public:
    Grid2D();
    int N;
    // underscore so we don't mix up input and class parameters
    Grid2D(int N_, int M_, double xmin_, double xmax_, double ymin_, double ymax_);
    // because dx and dy are private we need function
    double get_dx() ;
    double get_dy();
    double get_xmin() ;
    double get_ymin() ;
    double get_xmax() ;
    double get_ymax() ;
    int get_N() ;
    int get_M() ;
    // these functionss move between vector and 2d space
    int i_from_n(int n) ;
    int j_from_n(int n) ;
    int n_from_ij(int i, int j) ;
    // these function give us phyical location on the grid
    double x_from_n(int n) ;
    double y_from_n(int n) ;
    // these are derivatives 1st and 2nd order (also forward and backward)
    double dx_forward(std::vector <double> & function, int n);
    double dx_backward(std::vector <double> & function, int n);
    double dy_forward(std::vector <double> & function, int n);
    double dy_backward(std::vector <double> & function, int n);
    double dx_2(std::vector <double> & function, int n);
    double dy_2(std::vector <double> & function, int n);


    void initialize_VTK_file(std::string file_name);
    void print_VTK_Format( std::vector<double> &F, std::string data_name, std::string file_name );
};

#endif // GRID2D_H
