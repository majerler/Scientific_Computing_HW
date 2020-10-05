#include <iostream>
#include <grid2d.h>
using namespace std;



int main()
{
    cout << "Hello World!" << endl;

    int N = 2;
    int M = 3;
    double xmin = 0;
    double xmax = 2;
    double ymin = 0;
    double ymax = 2;

    // This command below creates an instance of the grid2d class
    Grid2D grid;

    //lets us fancy constuctor
    Grid2D new_grid(N, M, xmin, xmax, ymin, ymax);

    cout << "grid resolution of y: " << new_grid.get_dy()  << endl;

    cout << "x cord n = 4: " << new_grid.x_from_n(4) << endl;

    return 0;
}
