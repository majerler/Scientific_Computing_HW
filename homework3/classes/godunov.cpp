#include "godunov.h"

godunov::godunov()
{

}

godunov::godunov(Grid2D grid_,  velocity_X * velx_, velocity_Y * vely_, std::vector <double> solution_, std::vector <double> original_)
{
    grid = grid_;
    sol = solution_;
    og = original_;
    velx = velx_;
    vely = vely_;

}


void godunov::upwind(double dt)
{
    //#pragma omp parallel for
    for (int n = 0; n < grid.get_M() * grid.get_N(); ++n)
    {
        double x = grid.x_from_n(n);
        double y = grid.y_from_n(n);

//        double v_x = (*velx)(x,y);
//        double v_y = (*vely)(x,y);

        double dx1, dy1;
        double dx2, dy2;

        double dxf = grid.dx_forward( sol, n);
        double dxb = grid.dx_backward( sol, n);
        double dyf = grid.dy_forward( sol, n);
        double dyb = grid.dy_backward( sol, n);

        if ( og[n]*dxb <= 0 & og[n]*dxf <= 0 ){
            dx1 = dxf ;
            dx2 = minmod( grid.dx_2(sol, n), grid.dx_2(sol, n+1) );}

        if ( og[n]*dxb >= 0 & og[n]*dxf >= 0 ){
            dx1 = dxb ;
            dx2 = minmod( grid.dx_2(sol, n), grid.dx_2(sol, n-1) );}

        if ( og[n]*dxb <= 0 & og[n]*dxf >= 0 ){
            dx1 = 0.;
            dx2 = 0.;}

        if ( og[n]*dxb >= 0 & og[n]*dxf <= 0 ){
            if( abs(dxb) >= abs(dxf)){
                dx1 = dxb ;
                dx2 = minmod( grid.dx_2(sol, n), grid.dx_2(sol, n-1) );}
            else{
                dx1 = dxf ;
                dx2 = minmod( grid.dx_2(sol, n), grid.dx_2(sol, n+1) );}
        }


        if ( og[n]*dyb <= 0 & og[n]*dyf <= 0 ){
            dy1 = dyf ;
            dy2 = minmod( grid.dy_2(sol, n), grid.dx_2(sol, n+grid.get_N()) );}

        if ( og[n]*dyb >= 0 & og[n]*dyf >= 0 ){
            dy1 = dyb ;
            dy2 = minmod( grid.dy_2(sol, n), grid.dx_2(sol, n-grid.get_N()) );}

        if ( og[n]*dyb <= 0 & og[n]*dyf >= 0 ){
            dy1 = 0. ;
            dy2 = 0.;}

        if ( og[n]*dyb >= 0 & og[n]*dyf <= 0 ){
            if( abs(dyb) >= abs(dyf)){
                dy1 = dyb ;
                dy2 = minmod( grid.dy_2(sol, n), grid.dx_2(sol, n-grid.get_N()) );}
            else{
                dy1 = dyf ;
                dy2 = minmod( grid.dy_2(sol, n), grid.dx_2(sol, n+grid.get_N()) );}
        }

        double eps = 0.000001;

        double v_x = dx1 / sqrt( dx1*dx1 + dy1*dy1 + eps*eps );
        double v_y = dy1 / sqrt( dx1*dx1 + dy1*dy1 + eps*eps );

//        if (og[n] > 0){
//            sol[n] = sol[n] + dt * ( 1 - (v_x * dx1+ v_y * dy1)  );
//        }
//        else{
//            sol[n] = sol[n] - dt * ( 1 - (v_x * dx1 + v_y * dy1)  );
//        }

        if (og[n] > 0){
            sol[n] = sol[n] + dt * ( 1 - ( v_x*( dx1 + dx2*grid.get_dx()/2 )+ v_y*( dy1 + dy2*grid.get_dy()/2) )  ) ;
        }
        else{
            sol[n] = sol[n] - dt * ( 1 - ( v_x*( dx1 + dx2*grid.get_dx()/2 )+ v_y*( dy1 + dy2*grid.get_dy()/2) )  ) ;
        }



    }
}


std::vector <double> godunov::get_sol()
{
    return sol;
}
