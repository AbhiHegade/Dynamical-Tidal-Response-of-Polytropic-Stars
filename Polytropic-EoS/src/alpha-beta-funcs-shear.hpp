#ifndef _ALPHA_BETA_FUNCS_SHEAR_HPP_
#define _ALPHA_BETA_FUNCS_SHEAR_HPP_

#include <vector>
#include <functional>
#include <Eigen/Dense>

void get_alpha_beta_funcs_ders_and_shear_source(double l, double Omega, double xi,
double theta, 
double n, double brel, 
double gamma_0,
double lambda, double nu,
double feta, double feta_der, double feta_der_2,
double H00, double V0, double W0, double H00der, std::vector<double> &dydx);

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




#endif