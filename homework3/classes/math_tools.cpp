#include <../classes/math_tools.h>
//#include <../classes/grid2d.h>

#include <math.h>

double bilinear_interpolation(Grid2D & grid, std::vector <double> & func, double x, double y)
{


    // Get x grid points
    int i_min = (int) floor( (x - grid.get_xmin()) / grid.get_dx() );
        i_min = std::max(0, i_min);
    int i_max = (int) ceil ( (x - grid.get_xmin()) / grid.get_dx() );
        i_max = std::min(grid.get_N()-1, i_max);

        if (i_min == i_max){

            if (i_min == 0){
                i_max = i_min+1;
            }
            else{
                i_min = i_max-1;
            }
        }


    // Get y grid points
    int j_min = (int) floor( (y - grid.get_ymin()) / grid.get_dy() );
        j_min = std::max(0, j_min);
    int j_max = (int) ceil ( (y - grid.get_ymin()) / grid.get_dy() );
        j_max = std::min(grid.get_M()-1, j_max);

        if (j_min == j_max){

            if (j_min == 0){
                j_max = j_min+1;
            }
            else{
                j_min = j_max-1;
            }
        }


    int corn_00 = grid.n_from_ij(i_min, j_min);
    int corn_01 = grid.n_from_ij(i_min, j_max);
    int corn_10 = grid.n_from_ij(i_max, j_min);
    int corn_11 = grid.n_from_ij(i_max, j_max);

    double x_min = grid.x_from_n(corn_00);
    double y_min = grid.y_from_n(corn_00);

    double x_max = grid.x_from_n(corn_11);
    double y_max = grid.y_from_n(corn_11);



    double dx = grid.get_dx();
    double dy = grid.get_dy();


    // Lets interpolate
    return  1. / (dx * dy) * ( func[corn_00] * (x_max - x) * (y_max - y) +
                                  func[corn_01] * (x_max - x) * (y - y_min) +
                                  func[corn_10] * (x - x_min) * (y_max - y) +
                                  func[corn_11] * (x - x_min) * (y - y_min) );


}


double minmod(double l, double m){
    if (l*m < 0) return 0;
    else if(abs(l) < abs(m)) return l;
    else return m;
    }


double eno_interpolation(Grid2D & grid, std::vector <double> & func, double x, double y)
{


    // Get x grid points
    int i_min = (int) floor( (x - grid.get_xmin()) / grid.get_dx() );
        i_min = std::max(0, i_min);
    int i_max = (int) ceil ( (x - grid.get_xmin()) / grid.get_dx() );
        i_max = std::min(grid.get_N()-1, i_max);

        if (i_min == i_max){

            if (i_min == 0){
                i_max = i_min+1;
            }
            else{
                i_min = i_max-1;
            }
        }


    // Get y grid points
    int j_min = (int) floor( (y - grid.get_ymin()) / grid.get_dy() );
        j_min = std::max(0, j_min);
    int j_max = (int) ceil ( (y - grid.get_ymin()) / grid.get_dy() );
        j_max = std::min(grid.get_M()-1, j_max);

        if (j_min == j_max){

            if (j_min == 0){
                j_max = j_min+1;
            }
            else{
                j_min = j_max-1;
            }
        }


    int corn_00 = grid.n_from_ij(i_min, j_min);
    int corn_01 = grid.n_from_ij(i_min, j_max);
    int corn_10 = grid.n_from_ij(i_max, j_min);
    int corn_11 = grid.n_from_ij(i_max, j_max);

    double x_min = grid.x_from_n(corn_00);
    double y_min = grid.y_from_n(corn_00);

    double x_max = grid.x_from_n(corn_11);
    double y_max = grid.y_from_n(corn_11);

    double dx = grid.get_dx();
    double dy = grid.get_dy();

    double dxx_00 = grid.dx_2(func, corn_00);
    double dxx_01 = grid.dx_2(func, corn_01);
    double dxx_10 = grid.dx_2(func, corn_10);
    double dxx_11 = grid.dx_2(func, corn_11);

    double dxx_0 = minmod(dxx_00, dxx_01);
    double dxx_1 = minmod(dxx_10, dxx_11);

    double dxx  = minmod(dxx_0,  dxx_1);

    double dyy_00 = grid.dy_2(func, corn_00);
    double dyy_01 = grid.dy_2(func, corn_01);
    double dyy_10 = grid.dy_2(func, corn_10);
    double dyy_11 = grid.dy_2(func, corn_11);

    double dyy_0 = minmod(dyy_00, dyy_01);
    double dyy_1 = minmod(dyy_10, dyy_11);

    double dyy = minmod(dyy_0, dyy_1);



    // Lets interpolate
    return  ( 1. / (dx * dy) * ( func[corn_00] * (x_max - x) * (y_max - y)        +
                                 func[corn_01] * (x_max - x) * (y - y_min)        +
                                 func[corn_10] * (x - x_min) * (y_max - y)        +
                                 func[corn_11] * (x - x_min) * (y - y_min) ) )    -
                               ( dxx           * (x - x_min) * (x_max - x) / 2 ) -
                               ( dyy           * (y - y_min) * (y_max - y) / 2 );


}



