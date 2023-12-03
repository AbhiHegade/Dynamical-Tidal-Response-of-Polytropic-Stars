#ifndef _CONSERVATIVE_TIDES_HPP_
#define _CONSERVATIVE_TIDES_HPP_

#include <vector>
#include <string>
#include <functional>

#include "alpha-beta-funcs.hpp"
#include "local_expansions.hpp"

class CONSERVATIVE_TIDE{
public:
    CONSERVATIVE_TIDE(double l, double Omega, 
    double n, double brel,
    double nu0,
    double gamma_0,
    double mu1,
    double xi1,
    std::function <double(double)> thetafunc, 
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc);
    ~CONSERVATIVE_TIDE(void);
//----------------------------------------------
    double l; 
    double Omega; 
    double n; double brel;
    double nu0;
    double gamma_0;
    double mu1;
    double xi1;
    std::function <double(double)> thetafunc;
    std::function <double(double)> lambdafunc;
    std::function <double(double)> nufunc;

    std::vector<double> H0_ext;
    std::vector<double> H1_ext;

    RHS_PF_SOURCELESS rhs;
    SERIES_ORIGIN_SOURCELESS series_origin;
    SERIES_SURF series_surf;
//--------------------------------------------
    void solve_and_match( int iterate, 
    double dr, double rl, double rmatch, double rr, 
    vector<double> &icsl_val, vector<double> &icsr_val,
    double &norm, double &k2   );

    void solve_and_iterate(int steps, double dr, double rl, double rmatch, double rr, double &k2, double &norm, std::vector<double> &ssl, std::vector<double> &ssr);

    void solve_and_write(
    int N_write, double rl,  double rr,
    std::vector<double> icsl, std::vector<double> icsr,
    std::vector<double> &H0_v,
    std::vector<double> &W0_v,
    std::vector<double> &V_v,
    std::vector<double> &H1_v,
    std::vector<double> &rvals_v);

    void solve_and_iterate_and_write(int N_write,
    int steps, double dr, double rl, double rmatch, double rr,
    double &k2,
    std::vector<double> &H0_v,
    std::vector<double> &W0_v,
    std::vector<double> &V_v,
    std::vector<double> &H1_v,
    std::vector<double> &rvals_v,
    std::vector<double> &vals_origin,
    double &norm_adiabat);




private:
int get_H0_flag;
    //variable order
    //H0,W0,H1,W1
//-------------------------------------------------------------------------
void get_H0_and_H1_ext();
void fitting_coeffs_origin(double r, 
double H0_origin, double W0_origin, 
double &H_0, double &W_0 );
// //=========================================================================
// void series_surf_0(
//         double r,
//         const std::vector<double> yin, std::vector<double> &y );
//==========================================================================
// void step_0(
//      state_type &y, double &r, double &dr);
//==========================================================================

};



#endif
