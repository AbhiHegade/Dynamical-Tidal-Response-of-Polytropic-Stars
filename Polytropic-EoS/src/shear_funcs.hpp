#ifndef _ALPHA_BETA_FUNCS_SHEAR_HPP_
#define _ALPHA_BETA_FUNCS_SHEAR_HPP_

#include <vector>
#include <functional>
#include <Eigen/Dense>
#include "alpha-beta-funcs.hpp"
#include "local_expansions.hpp"

void get_alpha_beta_funcs_ders_and_shear_source(double l, double Omega, double xi,
double theta, 
double n, double brel, 
double gamma_0,
double lambda, double nu,
double feta, double feta_der, double feta_der_2,
double H00, double W0, double V0,  double H00der, 
std::vector<double> y,
std::vector<double> &dydx);

//======================================================================================================
class RHS_PF_SHEAR{
    
public:
    RHS_PF_SHEAR(
    double l, double Omega, 
    double n, double brel,
    double gamma_0,
    std::function <double(double)> thetafunc,  
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc,
    std::function <double(double, double)> fetafunc,
    std::function <double(double, double)> fetaderfunc,
    std::function <double(double, double)> fetader2func,
    std::function <double(double)> H00func,
    std::function <double(double)> W0func,
    std::function <double(double)> V0func,
    std::function <double(double)> H00derfunc);

    ~RHS_PF_SHEAR(void);
//----------------------------------------------
    double l;
    double Omega;
    double n;
    double brel;
    double gamma_0;
    std::function <double(double)> thetafunc;
    std::function <double(double)> lambdafunc; 
    std::function <double(double)> nufunc;
//---------------------------------------------
    std::function <double(double, double)> fetafunc;
    std::function <double(double, double)> fetaderfunc;
    std::function <double(double, double)> fetader2func;

    std::function <double(double)> H00func;
    std::function <double(double)> W0func;
    std::function <double(double)> V0func;
    std::function <double(double)> H00derfunc;
//---------------------------------------------
    void operator()( const std::vector<double> &y , std::vector<double> &dydt , double r );

private:
//============================================================================
    void get_dydx(double xi, const std::vector<double> &y, std::vector<double> &dydx);

    
};

//======================================================================================================
class SERIES_ORIGIN_SHEAR{

    public:
    SERIES_ORIGIN_SHEAR(double l, double Omega, double nu0, double n, double brel, double gamma_0, 
    double feta_c, double feta_der_c, double feta_der2_c, double H0_0, double W0_0);
    ~SERIES_ORIGIN_SHEAR(void);
//----------------------------------------------
    double l; 
    double Omega; 
    double nu0;
    double n; 
    double brel; 
    double gamma_0;
    double feta_c;
    double feta_der_c;
    double feta_der2_c;
    double H0_0;
    double W0_0;

    void operator()( double xi, const std::vector<double> yin, std::vector<double> &y );


    private:


};

//==================================================================================

class DISS_SHEAR_TIDE{
public:
    DISS_SHEAR_TIDE(double l, double Omega, 
    double n, double brel,
    double nu0,
    double gamma_0,
    double mu1,
    double xi1,
    std::function <double(double)> thetafunc, 
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc,
    std::function <double(double, double)> fetafunc,
    std::function <double(double, double)> fetaderfunc,
    std::function <double(double, double)> fetader2func,
    std::function <double(double)> H00, 
    std::function <double(double)> W0, 
    std::function <double(double)> V0,
    std::function <double(double)> H00der,
    double H0_0, double W0_0
        );
    ~DISS_SHEAR_TIDE(void);
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

    std::function <double(double, double)> fetafunc;
    std::function <double(double, double)> fetaderfunc;
    std::function <double(double, double)> fetader2func;
    std::function <double(double)> H00func;
    std::function <double(double)> W0func;
    std::function <double(double)> V0func;
    std::function <double(double)> H00derfunc;

    double H0_0; double W0_0;

    RHS_PF_SHEAR rhspart;
    SERIES_ORIGIN_SHEAR series_origin_part;
    RHS_PF_SOURCELESS rhshom;
    SERIES_ORIGIN_SOURCELESS series_origin_hom;
    SERIES_SURF series_surf;
//--------------------------------------------
    void solve_and_match( int iterate, 
    double dr, double rl, double rmatch, double rr, 
    vector<double> &icsl_val, vector<double> &icsr_val,
    double &norm, double &k2   );

    void solve_and_iterate(int steps, double dr, double rl, double rmatch, double rr, double &k2, double &norm);

private:
int get_H0_flag;
    //variable order
    //H0,W0,H1,W1
//-------------------------------------------------------------------------
void get_H0_and_H1_ext();


};


#endif