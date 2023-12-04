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
// #include <random>
// const double E = 2.718281828459045;
// const double Pi = 3.14159265358979323846;
// const double lower_bound = 0.5;
// const double upper_bound = 1.;
 
// std::uniform_real_distribution<double> unif_g(lower_bound,
//                                            upper_bound);
// std::default_random_engine re_g;
//-------------------------------------------------
#include "alpha-beta-funcs.hpp"
#include "local_expansions.hpp"
#include "conservative_tides.hpp"
#include "ode_steppers.hpp"
#include "outputfiles.hpp"
//-------------------------------------------------
CONSERVATIVE_TIDE::CONSERVATIVE_TIDE(double l, double Omega, 
    double n, double brel,
    double nu0,
    double gamma_0,
    double mu1,
    double xi1,
    std::function <double(double)> thetafunc, 
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc):
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
rhs(l, Omega, n, brel,gamma_0, thetafunc, lambdafunc, nufunc),
series_origin(l, Omega, nu0, n, brel, gamma_0),
series_surf(l, Omega, n, brel, gamma_0, mu1, xi1),
get_H0_flag{0}
{

}
CONSERVATIVE_TIDE::~CONSERVATIVE_TIDE(void)
{

}

//-------------------------------------------------------------------------------------------------------------------------------------
void CONSERVATIVE_TIDE::get_H0_and_H1_ext( )
{
    if(get_H0_flag ==0)
    {   
        obtain_H0_and_H1_ext(l, Omega, n, brel, mu1, xi1, H0_ext, H1_ext );

        get_H0_flag = 1;
    }
    else{

    }
    

}
//-------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------
void CONSERVATIVE_TIDE::solve_and_match(int iterate, 
double dr, double rl, double rmatch, double rr, 
vector<double> &icsl_val, vector<double> &icsr_val,
double &norm, double &k2 )
{
    double Rstar = xi1;
    assert(rl>0);
    assert(rr<Rstar);
    assert(rl<rmatch);
    assert(rr>rmatch);

    // Left solutions
    //-----------------------------------------
    int len_var = 4;

    vector<double> yl_hom1(len_var), yl_hom2(len_var);
    vector<double> ics1{0.1,0.}, ics2{0.1,0.1} ;
    vector<double> ssl_hom1(len_var), ssl_hom2(len_var);


    if(iterate)
    {
        ssl_hom1 = icsl_val;

        solve_adaptive(rhs, series_origin,  series_surf, "left", 0, rl, rmatch, Rstar, ics1, ssl_hom1, yl_hom1, dr);
        
    }
    else{
        solve_adaptive(rhs, series_origin,  series_surf, "left", 1, rl, rmatch, Rstar, ics1, ssl_hom1, yl_hom1, dr);
    }
    
    solve_adaptive(rhs, series_origin,  series_surf, "left", 1, rl, rmatch, Rstar, ics2, ssl_hom2, yl_hom2, dr);
  
    //-----------------------------------------
    // Right calculations

    vector<double> yr_hom1(len_var), yr_hom2(len_var), yr_hom3(len_var);
    vector<double> icsr1{0.1,0.,0.}, icsr2{0.1,0.1,0.}, icsr3{0.1,0.1,0.1} ;
    vector<double> ssr_hom1(len_var), ssr_hom2(len_var), ssr_hom3(len_var);

    
    if(iterate)
    {
        ssr_hom1 = icsr_val;
        solve_adaptive(rhs, series_origin,  series_surf, "right", 0, rr, rmatch, Rstar, ics1, ssr_hom1, yr_hom1, dr);
    }
    else{
          solve_adaptive(rhs, series_origin,  series_surf, "right", 1, rr, rmatch, Rstar, icsr1, ssr_hom1, yr_hom1, dr);
    }
    
    solve_adaptive(rhs, series_origin,  series_surf, "right", 1, rr, rmatch, Rstar, icsr2, ssr_hom2, yr_hom2, dr);
    solve_adaptive(rhs, series_origin,  series_surf, "right", 1, rr, rmatch, Rstar, icsr3, ssr_hom3, yr_hom3, dr);

    // //========================================================
    get_H0_and_H1_ext();

    Eigen::VectorXd b(6);
    b[0] = 0.; b[1] = 0.; b[2] = 0.; b[3] = 0.; 
    b[4] = H0_ext[0]; b[5] = H1_ext[0];

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
        icsl_val[i] = x[0]*ssl_hom1[i] + x[1]*ssl_hom2[i] ;
        
        icsr_val[i] = x[2]*ssr_hom1[i] + x[3]*ssr_hom2[i] + x[4]*ssr_hom3[i] ;

        
    }
   

    k2 = x[5];

    cout<<"k2_adiabat = "<<k2<<endl;

    solve_adaptive(rhs, series_origin,  series_surf, "left", 0, rl, rmatch, Rstar, ics1, icsl_val, yl_hom1, dr);
    solve_adaptive(rhs, series_origin,  series_surf, "right", 0, rr, rmatch, Rstar, ics1, icsr_val, yr_hom1, dr);

    norm = 0;
    for(int i = 0; i<len_var; ++i)
    {   
        norm += pow((yl_hom1[i] - yr_hom1[i]),2.);
    }
    norm = pow(norm,0.5);
    cout<<"norm = "<<norm<<endl;
    cout<<"--------------------------"<<endl;

}


//-------------------------------------------------

void CONSERVATIVE_TIDE::solve_and_iterate(int steps, double dr, double rl, double rmatch, double rr, double &k2, double &norm, std::vector<double> &ssl, std::vector<double> &ssr)
{

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
void CONSERVATIVE_TIDE::solve_and_write(
    int N_write, double rl, double rr,
    std::vector<double> icsl, std::vector<double> icsr,
    std::vector<double> &H0_v,
    std::vector<double> &W0_v,
    std::vector<double> &V_v,
    std::vector<double> &H1_v,
    std::vector<double> &rvals_v)
{
    RHS_PF_SOURCELESS rhs( l, Omega, n, brel, gamma_0, thetafunc,  lambdafunc, nufunc);
    double omega = Omega;

    assert(N_write %2 == 1);
    int N = (N_write+1)/2;

    double dx_new = (rr-rl)/(2*N - 2);
    // double drl = (rmatch-rl)/(N-1);
    // double drr = (rr-rmatch)/(N-1);
    double drl = dx_new;
    double drr = dx_new;

    double drl_start = drl/4.;
    double drr_start = drr/4.;
    double rl_start = rl;
    double rr_start = rr;

    cout<<"omega = "<<omega<<endl;
    vector<double> H0_l(N),H1_l(N),W0_l(N), V_l(N);
    vector<double> H0_r(N-1),H1_r(N-1),W0_r(N-1), V_r(N-1);
    vector<double> rl_v(N), rr_v(N-1);


    {
        int i=0;

        H0_l[i] = icsl[0];
        W0_l[i] = icsl[1];
        V_l[i] = icsl[2];
        H1_l[i] = icsl[3];

        H0_r[N-2-i] = icsr[0];
        W0_r[N-2-i] = icsr[1];
        V_r[N-2-i] = icsr[2];
        H1_r[N-2-i] = icsr[3];

        rl_v[i] = rl;
        rr_v[N-2-i] = rr;

    }
    for(int i =1; i<N-1; ++i)
    {   


        step_controlled(rhs, "left",icsl, rl_start, drl_start, drl); 

        step_controlled(rhs, "right", icsr, rr_start, drr_start, drr);

        rl_v[i] = rl_start;
        rr_v[N-2-i] = rr_start;

        H0_l[i] = icsl[0];
        W0_l[i] = icsl[1];
        V_l[i] = icsl[2];
        H1_l[i] = icsl[3];

        H0_r[N-2-i] = icsr[0];
        W0_r[N-2-i] = icsr[1];
        V_r[N-2-i] = icsr[2];
        H1_r[N-2-i] = icsr[3];
    }
    {
        int i = N-1;
        
        step_controlled(rhs, "left",icsl, rl_start, drl_start, drl); 

        step_controlled(rhs, "right", icsr, rr_start, drr_start, drr);

        // step_controlled_0("left",icsl, rl_start, drl_start, drl); 

        // step_controlled_0("right", icsr, rr_start, drr_start, drr);

        rl_v[i] = rl_start;
        
        H0_l[i] = icsl[0];
        W0_l[i] = icsl[1];
        V_l[i] = icsl[2];
        H1_l[i] = icsl[3];

        double res = 0.;
        for(int i=0; i<4; ++i)
        {
            res += pow(icsl[i]-icsr[i],2.);
        }
        cout<<"res_interp = "<<pow(res,0.5)<<endl;

    }

    combine_vec(H0_l,H0_r,H0_v);
    combine_vec(W0_l,W0_r,W0_v);
    combine_vec(V_l,V_r,V_v);
    combine_vec(H1_l,H1_r,H1_v);
    combine_vec(rl_v,rr_v,rvals_v);
}

//-------------------------------------------------

void CONSERVATIVE_TIDE::fitting_coeffs_origin(double r, 
double H0_origin, double W0_origin, 
double &H_0, double &W_0 )
{
    Eigen::MatrixXd A(2,2);
    Eigen::VectorXd b_v(2);
    b_v[0] = H0_origin; b_v[1] = W0_origin;

    double xi = r;

    get_matrix_ord_2_origin_just_H_W(l, Omega, nu0, n, brel, gamma_0, xi, A);
//----------------------------------
    Eigen::VectorXd x(2);
    x = A.fullPivLu().solve(b_v);

    H_0 = x[0];
    W_0 = x[1];

}

//-------------------------------------------------------------

void CONSERVATIVE_TIDE::solve_and_iterate_and_write(int N_write,
    int steps, double dr, double rl, double rmatch, double rr,
    double &k2,
    std::vector<double> &H0_v,
    std::vector<double> &W0_v,
    std::vector<double> &V_v,
    std::vector<double> &H1_v,
    std::vector<double> &rvals_v,
    std::vector<double> &vals_origin,
    double &norm_adiabat)
{

    int len_var = 4;
    // double norm = 0.;
    vector<double> ssl(len_var);
    vector<double> ssr(len_var);
    ssl[0] = 1;
    ssr[0] = 1;
    solve_and_iterate(steps, dr, rl, rmatch, rr, k2, norm_adiabat, ssl, ssr);

    solve_and_write(
    N_write,
    rl,  
    rr,
    ssl,ssr,
    H0_v,
    W0_v,
    V_v,
    H1_v,
    rvals_v);


    fitting_coeffs_origin(rl, 
    ssl[0], ssl[1], 
    vals_origin[0],vals_origin[1]);

    out_put_vec(ssl);
    cout<<"---"<<endl;
    out_put_vec(ssr);
    cout<<"---"<<endl;
    out_put_vec(vals_origin);
    cout<<"-------"<<endl;

    cout<<H0_v.back()<<endl;
    cout<<W0_v.back()<<endl;
    cout<<V_v.back()<<endl;
    cout<<H1_v.back()<<endl;
    
    
}
