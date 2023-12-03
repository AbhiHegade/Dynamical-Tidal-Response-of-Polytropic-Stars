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



#endif