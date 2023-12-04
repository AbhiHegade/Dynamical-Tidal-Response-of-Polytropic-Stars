#ifndef _ALPHA_BETA_FUNCS_HPP_
#define _ALPHA_BETA_FUNCS_HPP_

#include <vector>
#include <functional>
//----------------------
void alpha_funcs_source_less(double l, double Omega, double xi, double theta, 
 double n, double brel, double lambda, double nu, std::vector<double> &alpha_vals);

void beta_funcs_and_ders_source_less(double l, double Omega, double xi,
double theta, 
double n, double brel, double lambda, double nu, std::vector<double> &beta_vals, std::vector<double> &beta_vals_der);

void alpha_H2_W_V_vals_source_less(double l, double Omega, double xi,
double theta, 
double n, double brel, 
double gamma_0,
double lambda, double nu,
std::vector<double> &alpha_H2_vals, std::vector<double> &alpha_W_vals, std::vector<double> &alpha_V_vals);
//==========================================================================
//==========================================================================
class RHS_PF_SOURCELESS{
    
public:
    RHS_PF_SOURCELESS(
    double l, double Omega, 
    double n, double brel,
    double gamma_0,
    std::function <double(double)> thetafunc,  
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc);
    ~RHS_PF_SOURCELESS(void);
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
    void operator()( const std::vector<double> &y , std::vector<double> &dydt , double xi );

private:
//============================================================================
    void get_dydx(double xi, const std::vector<double> &y, std::vector<double> &dydx);

    
};

//==========================================================================
//-----------------Source and functions for bulk viscosity------------------------
void alpha_funcs_source_bulk(double l, double Omega, double xi, double theta, 
 double n, double brel, double lambda, double nu, double &alpha_4, double &alpha_7);

void beta_funcs_and_ders_source_bulk(double l, double Omega, double xi, double theta, 
 double n, double brel, double lambda, double nu, 
double &beta_4, double &beta_7,
double &beta_der_4, double &beta_der_7);

void alpha_H2_W_V_vals_source_bulk(double l, double Omega, double xi,
double theta, 
double n, double brel, 
double gamma_0,
double lambda, double nu,
double &alphaH_4, double &alphaH_7, double &alphaH_8, 
double &alphaW_4, double &alphaW_7,  
double &alphaV_4, double &alphaV_7, double &alphaV_8 );

//======================================================================================================
double S_Omega_bulk(double Omega, double xi,
double theta, 
double n, double brel, 
double gamma_0,
double lambda, double nu,
double fzeta,
double H0,
double V,
double W
);
//======================================================================================================
class RHS_PF_BULK{
    
public:
    RHS_PF_BULK(
    double l, double Omega, 
    double n, double brel,
    double gamma_0,
    std::function <double(double)> thetafunc,  
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc,
    std::function <double(double, double)> fzetafunc,
    std::function <double(double)> H00func,
    std::function <double(double)> W0func,
    std::function <double(double)> V0func);

    ~RHS_PF_BULK(void);
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
    std::function <double(double, double)> fzetafunc;

    std::function <double(double)> H00func;
    std::function <double(double)> W0func;
    std::function <double(double)> V0func;
//---------------------------------------------
    void operator()( const std::vector<double> &y , std::vector<double> &dydt , double r );

private:
//============================================================================
    void get_dydx(double xi, const std::vector<double> &y, std::vector<double> &dydx);

    
};

//======================================================================================================











#endif