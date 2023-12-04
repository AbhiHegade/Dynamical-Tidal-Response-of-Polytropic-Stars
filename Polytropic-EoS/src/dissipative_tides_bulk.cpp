#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <fstream>
#include <Eigen/Dense>
using Eigen::seq;
#include <functional>
#include <cmath>
using std::pow;
using std::abs;
using std::log;
#include<cassert>
//--------------------------------------------------
#include "alpha-beta-funcs.hpp"
#include "local_expansions.hpp"
#include "dissipative_tides_bulk.hpp"
#include "ode_steppers.hpp"
#include "outputfiles.hpp"
//-------------------------------------------------
DISS_BULK_TIDE::DISS_BULK_TIDE(double l, double Omega, 
    double n, double brel,
    double nu0,
    double gamma_0,
    double mu1,
    double xi1,
    std::function <double(double)> thetafunc, 
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc,
    std::function <double(double, double)> fzetafunc,
    std::function <double(double)> H00, 
    std::function <double(double)> W0, 
    std::function <double(double)> V0,
    double H0_0, double W0_0):
l{l},
Omega{Omega},
n{n},
brel{brel},
nu0{nu0},
gamma_0{gamma_0},
mu1{mu1},
xi1{xi1},
thetafunc{thetafunc},
lambdafunc{lambdafunc},
nufunc{nufunc},
H0_ext(2),
H1_ext(2),
fzetafunc{fzetafunc},
H00func{H00},
W0func{W0},
V0func{V0},
H0_0{H0_0},
W0_0{W0_0},
rhspart(l, Omega, n, brel, gamma_0, thetafunc, lambdafunc, nufunc, fzetafunc, H00func, W0func, V0func),
series_origin_part(l, Omega, nu0, n, brel, gamma_0, fzetafunc(1.0,n), H0_0, W0_0),
rhshom(l, Omega, n, brel,gamma_0, thetafunc, lambdafunc, nufunc),
series_origin_hom(l, Omega, nu0, n, brel, gamma_0),
series_surf(l, Omega, n, brel, gamma_0, mu1, xi1),
get_H0_flag{0}
{
    
}
DISS_BULK_TIDE::~DISS_BULK_TIDE(void)
{

}

//-------------------------------------------------------------------------------------------------------------------------------------
void DISS_BULK_TIDE::get_H0_and_H1_ext( )
{
    if(get_H0_flag ==0)
    {
        obtain_H0_and_H1_ext(l, Omega, n, brel, mu1, xi1, H0_ext, H1_ext );
        
        H0_ext[0] = 0.;

        H1_ext[0] = 0.;

        get_H0_flag = 1;
    }
    else{

    }
    

}
//-------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------
void DISS_BULK_TIDE::solve_and_match(int iterate, 
double dr, double rl, double rmatch, double rr, 
vector<double> &icsl_val, vector<double> &icsr_val,
double &norm, double &k2 )
{   
    double omega = Omega;
    double Rstar = xi1;
    assert(rl>0);
    assert(rr<Rstar);
    assert(rl<rmatch);
    assert(rr>rmatch);

    // Left solutions
    //-----------------------------------------
    int len_var = 4;

    vector<double> yl_part(len_var), yl_hom1(len_var), yl_hom2(len_var);
    vector<double> ics1{1.,0.}, ics2{1.,1.} ;
    vector<double> ssl_part(len_var), ssl_hom1(len_var), ssl_hom2(len_var);


    if(iterate)
    {
        ssl_part = icsl_val;

        solve_adaptive(rhspart, series_origin_part,  series_surf, "left", 0, rl, rmatch, Rstar, ics1, ssl_part, yl_part, dr);
        
    }
    else{
         solve_adaptive(rhspart, series_origin_part,  series_surf, "left", 1, rl, rmatch, Rstar, ics1, ssl_part, yl_part, dr);
    }
    
    solve_adaptive(rhshom, series_origin_hom,  series_surf, "left", 1, rl, rmatch, Rstar, ics1, ssl_hom1, yl_hom1, dr);
    solve_adaptive(rhshom, series_origin_hom,  series_surf, "left", 1, rl, rmatch, Rstar, ics2, ssl_hom2, yl_hom2, dr);
  
    //-----------------------------------------
    // Right calculations

    vector<double> yr_part(len_var), yr_hom1(len_var), yr_hom2(len_var), yr_hom3(len_var);
    vector<double> icsr1{0.1,0.,0.}, icsr2{0.1,0.1,0.}, icsr3{0.1,0.1,0.1} ;
    vector<double> ssr_part(len_var), ssr_hom1(len_var), ssr_hom2(len_var), ssr_hom3(len_var);
    
    if(iterate)
    {
        ssr_part = icsr_val;
        solve_adaptive(rhspart, series_origin_part,  series_surf, "right", 0, rr, rmatch, Rstar, icsr1, ssr_part, yr_part, dr);
    }
    else{
          solve_adaptive(rhspart, series_origin_part,  series_surf, "right", 1, rr, rmatch, Rstar, icsr1, ssr_hom1, yr_hom1, dr);
    }
    
    solve_adaptive(rhshom, series_origin_hom,  series_surf, "right", 1, rr, rmatch, Rstar, icsr1, ssr_hom1, yr_hom1, dr);
    solve_adaptive(rhshom, series_origin_hom,  series_surf, "right", 1, rr, rmatch, Rstar, icsr2, ssr_hom2, yr_hom2, dr);
    solve_adaptive(rhshom, series_origin_hom,  series_surf, "right", 1, rr, rmatch, Rstar, icsr3, ssr_hom3, yr_hom3, dr);

    // //========================================================
    get_H0_and_H1_ext();

    Eigen::VectorXd b(6);

    for(int i=0; i<4; ++i)
    {
        b[i] = yr_part[i] - yl_part[i]; 
    }
    b[4] = H0_ext[0] - ssr_part[0]; 
    b[5] = H1_ext[0] - ssr_part[3];

    Eigen::MatrixXd A(6,6);

    for(int i=0; i<len_var; ++i)
    {
        A(i,0) = yl_hom1[i]; 
        A(i,1) = yl_hom2[i];

        A(i,2) = -yr_hom1[i]; 
        A(i,3) = -yr_hom2[i];
        A(i,4) = -yr_hom3[i];

    }
    A(4,2) = ssr_hom1[0]; A(4,3) = ssr_hom2[0]; 
    A(4,4) = ssr_hom3[0]; A(4,5) = -H0_ext[1];

    A(5,2) = ssr_hom1[3]; 
    A(5,3) = ssr_hom2[3];
    A(5,4) = ssr_hom3[3];
    A(5,5) = -H1_ext[1];
    


    Eigen::VectorXd x(6);
    x = A.fullPivLu().solve(b);
    vector<double> Cs(x.data(), x.data() + x.rows() * x.cols());
    check_field_isfinite(rmatch, Cs, "solution for matching coefficients in nan");
    double relative_error = (A*x - b).norm() / b.norm(); // norm() is L2 norm
    std::cout << "The relative error is:\n" << relative_error << std::endl;
    std::cout << "The determinant of A is " << A.determinant() << std::endl;

    
    for(int i =0; i<len_var; ++i)
    {   
        icsl_val[i] = ssl_part[i] + x[0]*ssl_hom1[i] + x[1]*ssl_hom2[i] ;
        
        icsr_val[i] = ssr_part[i] + x[2]*ssr_hom1[i] + x[3]*ssr_hom2[i] + x[4]*ssr_hom3[i] ;

    }
   

    k2 = x[5]/omega;

    cout<<"k2_tau = "<<k2<<endl;

    solve_adaptive(rhspart, series_origin_part,  series_surf, "left", 0, rl, rmatch, Rstar, ics1, icsl_val, yl_part, dr);
    solve_adaptive(rhspart, series_origin_part,  series_surf, "right", 0, rr, rmatch, Rstar, ics1, icsr_val, yr_part, dr);

    norm = 0;
    for(int i = 0; i<len_var; ++i)
    {   
        norm += pow((yl_part[i] - yr_part[i]),2.);
    }
    norm = pow(norm,0.5);
    cout<<"norm = "<<norm<<endl;
    cout<<"--------------------------"<<endl;

}


//-------------------------------------------------

void DISS_BULK_TIDE::solve_and_iterate(int steps, double dr, double rl, double rmatch, double rr, double &k2, double &norm)
{
    int len_var = 4;
    vector<double> ssl(len_var);
    vector<double> ssr(len_var);
    norm = 0.;
    ssl[0] = 1;
    ssr[0] = 1;
    solve_and_match(0, dr, rl, rmatch, rr, ssl, ssr, norm, k2);
    cout<<"norm = "<<norm<<endl;
    
    int counter = 1;
    while(counter<steps && norm > 1e-10)
    {   
        
        cout<<"===================================================="<<endl;
        cout<<"counter = "<<counter<<endl;
        cout<<"norm = "<<norm<<endl;
        solve_and_match(1, dr, rl, rmatch, rr, ssl, ssr, norm, k2);

        counter +=1;
    }
}
//-------------------------------------------------
