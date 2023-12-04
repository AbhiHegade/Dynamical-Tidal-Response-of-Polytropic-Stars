#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;
#include <cmath>
using std::pow;
using std::abs;
using std::log;
#include<cassert>
// #include <Eigen/Dense>
const double E = 2.718281828459045;
const double Pi = 3.14159265358979323846;
//-------------------------------------------------
#include "alpha-beta-funcs.hpp"
//-------------------------------------------------
void alpha_funcs_source_less(double l, double Omega, double xi, double theta, 
 double n, double brel, double lambda, double nu, std::vector<double> &alpha_vals)
{
double nu_der_1 = (-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*pow(xi,2)*pow(theta,1 + n) + 2*pow(brel,2)*pow(E,lambda)*n*pow(xi,2)*pow(theta,1 + n))/xi;

double denom = 
pow(E,nu)*(pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3)) + 4*brel*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*pow(E,lambda)*(-1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*(1 + n)*nu_der_1*pow(xi,3)*pow(Omega,2) + brel*(1 + n)*pow(nu_der_1,2)*pow(xi,4)*pow(Omega,2) + 4*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*pow(xi,4)*pow(Omega,4))
;
//----------------------
double alpha_val_0 = -2*pow(E,nu)*(pow(E,2*lambda + nu)*l*(1 + l) - pow(E,nu)*l*pow(1 + l,2) + pow(E,lambda + nu)*l*(-2 + 3*pow(l,2) + pow(l,3)) - 2*brel*(-7 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) - 2*brel*pow(E,lambda)*(3 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) - 2*brel*pow(E,lambda)*(1 + n)*pow(xi,2)*(pow(E,nu)*l*(1 + l) - 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,n) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*(2*pow(E,lambda + nu)*l*(1 + l) + pow(E,nu)*l*(-3 - 2*l + pow(l,2)) - 4*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n) + 4*pow(brel,4)*pow(E,2*lambda + nu)*l*(1 + l)*pow(1 + n,2)*pow(xi,4)*pow(theta,2 + 2*n));
alpha_vals[0]=alpha_val_0/denom ;
//----------------------
double alpha_val_1 = 4*brel*pow(E,lambda/2. + nu)*pow(1 + n,2)*pow(xi,2)*pow(Omega,2)*pow(theta,n)*(1 + brel*theta)*(-3 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(theta,1 + n));
alpha_vals[1]=alpha_val_1/denom ;
//----------------------
double alpha_val_2 = 4*brel*pow(E,lambda)*pow(1 + n,2)*pow(xi,2)*pow(Omega,2)*(-(pow(E,nu)*l*(1 + l)) + 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,n)*(1 + brel*theta);
alpha_vals[2]=alpha_val_2/denom ;
//----------------------
double alpha_val_3 = -2*pow(E,nu)*xi*(-(pow(E,nu)*l*(1 + l)) + pow(E,lambda + nu)*l*(1 + l) - 4*brel*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda + nu)*l*(1 + l)*(1 + n)*pow(xi,2)*pow(theta,1 + n));
alpha_vals[3]=alpha_val_3/denom ;
//----------------------
}
//=========================================================================================
void beta_funcs_and_ders_source_less(double l, double Omega, double xi,
double theta,  
double n, double brel, double lambda, double nu, std::vector<double> &beta_vals, std::vector<double> &beta_vals_der)
{

//-----------------------------------------------
double theta_der_1 = -0.5*((1 + brel*theta)*(-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(theta,1 + n)))/(brel*(1 + n)*xi);

double lambda_der_1 = (1 - pow(E,lambda) + 2*brel*pow(E,lambda)*pow(xi,2)*pow(theta,n) + 2*brel*pow(E,lambda)*n*pow(xi,2)*pow(theta,n))/xi;

double nu_der_1 = (-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*pow(xi,2)*pow(theta,1 + n) + 2*pow(brel,2)*pow(E,lambda)*n*pow(xi,2)*pow(theta,1 + n))/xi;


double nu_der_2 = pow(xi,-2) - pow(E,2*lambda)/pow(xi,2) + brel*pow(E,lambda)*(1 + pow(E,lambda))*(1 + n)*pow(theta,n) - pow(brel,2)*pow(E,lambda)*(-5 + 3*pow(E,lambda))*(1 + n)*pow(theta,1 + n) + 2*pow(brel,3)*pow(E,2*lambda)*pow(1 + n,2)*pow(xi,2)*pow(theta,1 + 2*n) - 2*pow(brel,4)*pow(E,2*lambda)*pow(1 + n,2)*pow(xi,2)*pow(theta,2 + 2*n)
; 
//-----------------------------
double denom = pow(E,nu)*xi*(pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3)) + 4*brel*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*pow(E,lambda)*(-1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*(1 + n)*nu_der_1*pow(xi,3)*pow(Omega,2) + brel*(1 + n)*pow(nu_der_1,2)*pow(xi,4)*pow(Omega,2) + 4*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*pow(xi,4)*pow(Omega,4));

double denom_der = pow(E,nu)*(pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3)) + pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3))*nu_der_1*xi + pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3))*(lambda_der_1 + nu_der_1)*xi + 12*brel*(1 + n)*pow(xi,2)*pow(Omega,2) - 12*brel*pow(E,lambda)*(-1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*pow(E,lambda)*(-1 + l + pow(l,2))*lambda_der_1*(1 + n)*pow(xi,3)*pow(Omega,2) - 12*brel*(1 + n)*nu_der_1*pow(xi,3)*pow(Omega,2) - 4*brel*pow(E,lambda)*(-1 + l + pow(l,2))*(1 + n)*nu_der_1*pow(xi,3)*pow(Omega,2) + brel*(1 + n)*pow(nu_der_1,2)*pow(xi,4)*pow(Omega,2) - 4*brel*(1 + n)*nu_der_2*pow(xi,4)*pow(Omega,2) + brel*(1 + n)*pow(nu_der_1,3)*pow(xi,5)*pow(Omega,2) + 2*brel*(1 + n)*nu_der_1*nu_der_2*pow(xi,5)*pow(Omega,2) + 20*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*pow(xi,4)*pow(Omega,4) + 4*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*(lambda_der_1 - nu_der_1)*pow(xi,5)*pow(Omega,4) + 4*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*nu_der_1*pow(xi,5)*pow(Omega,4))
;
//----------------------
double beta_val_0 = -2*(pow(E,2*nu)*pow(l,2)*pow(1 + l,2) + pow(E,lambda + 2*nu)*l*(2 + l - 3*pow(l,2) - 2*pow(l,3)) + pow(E,2*(lambda + nu))*l*(-2 - 2*l + pow(l,2) + pow(l,3)) + 4*brel*pow(E,nu)*(-3 - 2*l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*pow(E,2*lambda + nu)*(2 + 3*l + 3*pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) + brel*pow(E,lambda + nu)*(14 + 7*l + 3*pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) + 4*pow(brel,2)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,4) + 2*brel*pow(E,lambda + nu)*(1 + n)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,n) - 2*pow(brel,2)*pow(E,lambda + nu)*(1 + n)*pow(xi,2)*(-(pow(E,lambda + nu)*(-2 + l)*l*pow(1 + l,2)) + pow(E,nu)*pow(l,2)*(-3 - 2*l + pow(l,2)) - brel*(7 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 3*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n) + 4*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,3)*pow(xi,6)*pow(Omega,2)*pow(theta,1 + 2*n) - 4*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,2)*pow(xi,4)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,2 + 2*n));
beta_vals[0]=beta_val_0/denom ;
//----------------------
double beta_val_1 = -4*brel*pow(E,lambda/2.)*pow(1 + n,2)*pow(xi,2)*pow(Omega,2)*pow(theta,n)*(1 + brel*theta)*(-3*pow(E,nu)*l - pow(E,lambda + nu)*(-2 + pow(l,2)) + 2*brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda + nu)*l*(1 + n)*pow(xi,2)*pow(theta,1 + n));
beta_vals[1]=beta_val_1/denom ;
//----------------------
double beta_val_2 = 4*brel*pow(E,lambda)*pow(1 + n,2)*pow(xi,2)*pow(Omega,2)*pow(theta,n)*(1 + brel*theta)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n));
beta_vals[2]=beta_val_2/denom ;
//----------------------
double beta_val_3 = 2*pow(E,nu)*xi*(-(pow(E,nu)*pow(l,2)*(1 + l)) - pow(E,lambda + nu)*l*(-2 - 2*l + pow(l,2) + pow(l,3)) - 2*brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l) + 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n));
beta_vals[3]=beta_val_3/denom ;
//----------------------
//----------------------
double beta_val_der_0 = -2*(2*pow(E,2*nu)*pow(l,2)*pow(1 + l,2)*nu_der_1 + 2*pow(E,2*(lambda + nu))*l*(-2 - 2*l + pow(l,2) + pow(l,3))*(lambda_der_1 + nu_der_1) + pow(E,lambda + 2*nu)*l*(2 + l - 3*pow(l,2) - 2*pow(l,3))*(lambda_der_1 + 2*nu_der_1) + 8*brel*pow(E,nu)*(-3 - 2*l + pow(l,2))*(1 + n)*xi*pow(Omega,2) - 2*brel*pow(E,2*lambda + nu)*(2 + 3*l + 3*pow(l,2))*(1 + n)*xi*pow(Omega,2) + 2*brel*pow(E,lambda + nu)*(14 + 7*l + 3*pow(l,2))*(1 + n)*xi*pow(Omega,2) + 4*brel*pow(E,nu)*(-3 - 2*l + pow(l,2))*(1 + n)*nu_der_1*pow(xi,2)*pow(Omega,2) + brel*pow(E,lambda + nu)*(14 + 7*l + 3*pow(l,2))*(1 + n)*(lambda_der_1 + nu_der_1)*pow(xi,2)*pow(Omega,2) - brel*pow(E,2*lambda + nu)*(2 + 3*l + 3*pow(l,2))*(1 + n)*(2*lambda_der_1 + nu_der_1)*pow(xi,2)*pow(Omega,2) + 16*pow(brel,2)*pow(E,lambda)*pow(1 + n,2)*pow(xi,3)*pow(Omega,4) + 4*pow(brel,2)*pow(E,lambda)*lambda_der_1*pow(1 + n,2)*pow(xi,4)*pow(Omega,4) + 2*brel*pow(E,lambda + nu)*n*(1 + n)*theta_der_1*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,-1 + n) + 4*brel*pow(E,lambda + nu)*(1 + n)*xi*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,n) + 2*brel*pow(E,lambda + nu)*(1 + n)*(lambda_der_1 + nu_der_1)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,n) - 2*pow(brel,2)*pow(E,lambda + nu)*pow(1 + n,2)*theta_der_1*pow(xi,2)*(-(pow(E,lambda + nu)*(-2 + l)*l*pow(1 + l,2)) + pow(E,nu)*pow(l,2)*(-3 - 2*l + pow(l,2)) - brel*(7 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 3*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,n) + 2*brel*pow(E,lambda + nu)*(1 + n)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l)*nu_der_1 + 2*brel*(-3 + pow(E,lambda) - 2*l)*(1 + n)*xi*pow(Omega,2) + brel*pow(E,lambda)*lambda_der_1*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,n) + 4*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,3)*(1 + 2*n)*theta_der_1*pow(xi,6)*pow(Omega,2)*pow(theta,2*n) - 4*pow(brel,2)*pow(E,lambda + nu)*(1 + n)*xi*(-(pow(E,lambda + nu)*(-2 + l)*l*pow(1 + l,2)) + pow(E,nu)*pow(l,2)*(-3 - 2*l + pow(l,2)) - brel*(7 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 3*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n) - 2*pow(brel,2)*pow(E,lambda + nu)*(1 + n)*(lambda_der_1 + nu_der_1)*pow(xi,2)*(-(pow(E,lambda + nu)*(-2 + l)*l*pow(1 + l,2)) + pow(E,nu)*pow(l,2)*(-3 - 2*l + pow(l,2)) - brel*(7 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 3*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n) - 2*pow(brel,2)*pow(E,lambda + nu)*(1 + n)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(-3 - 2*l + pow(l,2))*nu_der_1 - pow(E,lambda + nu)*(-2 + l)*l*pow(1 + l,2)*(lambda_der_1 + nu_der_1) - 2*brel*(7 + 2*l)*(1 + n)*xi*pow(Omega,2) + 6*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*xi*pow(Omega,2) + 3*brel*pow(E,lambda)*(1 + l + pow(l,2))*lambda_der_1*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n) + 24*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,3)*pow(xi,5)*pow(Omega,2)*pow(theta,1 + 2*n) + 4*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,3)*(2*lambda_der_1 + nu_der_1)*pow(xi,6)*pow(Omega,2)*pow(theta,1 + 2*n) - 8*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,3)*theta_der_1*pow(xi,4)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + 2*n) - 4*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,2)*pow(xi,4)*(pow(E,nu)*pow(l,2)*(1 + l)*nu_der_1 + 2*brel*(1 + n)*xi*pow(Omega,2))*pow(theta,2 + 2*n) - 16*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,2)*pow(xi,3)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,2 + 2*n) - 4*pow(brel,4)*pow(E,2*lambda + nu)*pow(1 + n,2)*(2*lambda_der_1 + nu_der_1)*pow(xi,4)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,2 + 2*n));
beta_vals_der[0]=beta_val_der_0/denom - (denom_der/(denom*denom))*beta_val_0 ; 
//----------------------
double beta_val_der_1 = 2*brel*pow(E,lambda/2.)*pow(1 + n,2)*xi*pow(Omega,2)*pow(theta,n)*(-2*brel*theta_der_1*xi*(-3*pow(E,nu)*l - pow(E,lambda + nu)*(-2 + pow(l,2)) + 2*brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda + nu)*l*(1 + n)*pow(xi,2)*pow(theta,1 + n)) - 4*(1 + brel*theta)*(-3*pow(E,nu)*l - pow(E,lambda + nu)*(-2 + pow(l,2)) + 2*brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda + nu)*l*(1 + n)*pow(xi,2)*pow(theta,1 + n)) - lambda_der_1*xi*(1 + brel*theta)*(-3*pow(E,nu)*l - pow(E,lambda + nu)*(-2 + pow(l,2)) + 2*brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda + nu)*l*(1 + n)*pow(xi,2)*pow(theta,1 + n)) - (2*n*theta_der_1*xi*(1 + brel*theta)*(-3*pow(E,nu)*l - pow(E,lambda + nu)*(-2 + pow(l,2)) + 2*brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda + nu)*l*(1 + n)*pow(xi,2)*pow(theta,1 + n)))/theta - 2*xi*(1 + brel*theta)*(-3*pow(E,nu)*l*nu_der_1 - pow(E,lambda + nu)*(-2 + pow(l,2))*(lambda_der_1 + nu_der_1) + 4*brel*pow(E,lambda)*(1 + n)*xi*pow(Omega,2) + 2*brel*pow(E,lambda)*lambda_der_1*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda + nu)*l*pow(1 + n,2)*theta_der_1*pow(xi,2)*pow(theta,n) + 4*pow(brel,2)*pow(E,lambda + nu)*l*(1 + n)*xi*pow(theta,1 + n) + 2*pow(brel,2)*pow(E,lambda + nu)*l*(1 + n)*(lambda_der_1 + nu_der_1)*pow(xi,2)*pow(theta,1 + n)));
beta_vals_der[1]=beta_val_der_1/denom - (denom_der/(denom*denom))*beta_val_1 ; 
//----------------------
double beta_val_der_2 = 4*brel*pow(E,lambda)*pow(1 + n,2)*xi*pow(Omega,2)*pow(theta,n)*(brel*theta_der_1*xi*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n)) + 2*(1 + brel*theta)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n)) + lambda_der_1*xi*(1 + brel*theta)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n)) + (n*theta_der_1*xi*(1 + brel*theta)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n)))/theta + xi*(1 + brel*theta)*(pow(E,nu)*pow(l,2)*(1 + l)*nu_der_1 + 2*brel*pow(E,lambda)*(1 + n)*xi*pow(Omega,2) - 2*brel*(3 + 2*l)*(1 + n)*xi*pow(Omega,2) + brel*pow(E,lambda)*lambda_der_1*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,3)*theta_der_1*pow(xi,4)*pow(Omega,2)*pow(theta,n) + 8*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,3)*pow(Omega,2)*pow(theta,1 + n) + 2*pow(brel,3)*pow(E,lambda)*lambda_der_1*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n)));
beta_vals_der[2]=beta_val_der_2/denom - (denom_der/(denom*denom))*beta_val_2 ; 
//----------------------
double beta_val_der_3 = 2*pow(E,nu)*(-(pow(E,nu)*pow(l,2)*(1 + l)) - pow(E,lambda + nu)*l*(-2 - 2*l + pow(l,2) + pow(l,3)) - 2*brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l) + 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n) + nu_der_1*xi*(-(pow(E,nu)*pow(l,2)*(1 + l)) - pow(E,lambda + nu)*l*(-2 - 2*l + pow(l,2) + pow(l,3)) - 2*brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l) + 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n)) + xi*(-(pow(E,nu)*pow(l,2)*(1 + l)*nu_der_1) - pow(E,lambda + nu)*l*(-2 - 2*l + pow(l,2) + pow(l,3))*(lambda_der_1 + nu_der_1) - 4*brel*(3 + 2*l)*(1 + n)*xi*pow(Omega,2) + 4*brel*pow(E,lambda)*(1 + l + pow(l,2))*(1 + n)*xi*pow(Omega,2) + 2*brel*pow(E,lambda)*(1 + l + pow(l,2))*lambda_der_1*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,2)*pow(E,lambda)*pow(1 + n,2)*theta_der_1*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l) + 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,n) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l)*nu_der_1 + 4*brel*(1 + n)*xi*pow(Omega,2))*pow(theta,1 + n) + 4*pow(brel,2)*pow(E,lambda)*(1 + n)*xi*(pow(E,nu)*pow(l,2)*(1 + l) + 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n) + 2*pow(brel,2)*pow(E,lambda)*lambda_der_1*(1 + n)*pow(xi,2)*(pow(E,nu)*pow(l,2)*(1 + l) + 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2))*pow(theta,1 + n)));
beta_vals_der[3]=beta_val_der_3/denom - (denom_der/(denom*denom))*beta_val_3 ; 
//----------------------
}
//=========================================================================================
void alpha_H2_W_V_vals_source_less(double l, double Omega, double xi,
double theta, 
double n, double brel, 
double gamma_0,
double lambda, double nu,
std::vector<double> &alpha_H2_vals, std::vector<double> &alpha_W_vals, std::vector<double> &alpha_V_vals)
{

double theta_der_1 = -0.5*((1 + brel*theta)*(-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(theta,1 + n)))/(brel*(1 + n)*xi);

double lambda_der_1 = (1 - pow(E,lambda) + 2*brel*pow(E,lambda)*pow(xi,2)*pow(theta,n) + 2*brel*pow(E,lambda)*n*pow(xi,2)*pow(theta,n))/xi;

double nu_der_1 = (-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*pow(xi,2)*pow(theta,1 + n) + 2*pow(brel,2)*pow(E,lambda)*n*pow(xi,2)*pow(theta,1 + n))/xi;
//-----------------------------------------------------------------------------------------------------------
vector<double> alpha_vals(4), beta_vals(4), beta_vals_der(4);

alpha_funcs_source_less(l, Omega, xi, theta, n, brel, lambda, nu, alpha_vals);
beta_funcs_and_ders_source_less(l, Omega,  xi, theta, n, brel,lambda, nu, beta_vals, beta_vals_der);

double alpha_0 = alpha_vals[0], alpha_1 = alpha_vals[1], alpha_2 = alpha_vals[2], alpha_3 = alpha_vals[3];
double beta_0 = beta_vals[0], beta_1 = beta_vals[1], beta_2 = beta_vals[2], beta_3 = beta_vals[3];
double beta_der_0 = beta_vals_der[0], beta_der_1 = beta_vals_der[1], beta_der_2 = beta_vals_der[2], beta_der_3 = beta_vals_der[3];

//----------------------
double alphaW_0 = -(pow(E,lambda/2.)*xi*(brel - brel*alpha_0 + 1/(gamma_0*theta)));
//----------------------
double alphaW_1 = (-1 + pow(E,lambda) + 2*brel*theta*(-((1 + l)*gamma_0) + brel*pow(E,lambda/2.)*gamma_0*pow(xi,2)*alpha_1 + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(theta,n)))/(2.*brel*gamma_0*xi*theta);
//----------------------
double alphaW_2 = -(pow(E,lambda/2.)*xi*((l*(1 + l))/pow(xi,2) - brel*alpha_2 - ((1 + n)*pow(Omega,2))/(pow(E,nu)*gamma_0*theta)));
//----------------------
double alphaW_3 = brel*pow(E,lambda/2.)*xi*alpha_3;
//----------------------

//----------------------
double alphaV_0 = -0.5*(pow(E,nu)*((-1 + pow(E,lambda))*(-1 + n*(-1 + gamma_0)) + 2*brel*(1 + n)*theta*(-(gamma_0*(2*(l + nu_der_1*xi) + l*alpha_0 + xi*beta_0)) + brel*pow(E,lambda)*(-1 + n*(-1 + gamma_0))*pow(xi,2)*pow(theta,n))))/(brel*pow(1 + n,2)*gamma_0*xi*pow(Omega,2)*theta);
//----------------------
double alphaV_1 = -(pow(E,lambda/2.)/xi) + (pow(E,nu)*((4*l*(1 + n)*alpha_1)/xi + 4*(1 + n)*beta_1 + ((-1 + n*(-1 + gamma_0))*pow(-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(theta,1 + n),2))/(pow(brel,2)*pow(E,lambda/2.)*gamma_0*pow(xi,3)*theta)))/(4.*pow(1 + n,2)*pow(Omega,2));
//----------------------
double alphaV_2 = (-1 + pow(E,lambda) - l + (pow(E,nu)*(l*alpha_2 + xi*beta_2))/((1 + n)*pow(Omega,2)))/xi + ((-1 + pow(E,lambda))*(-1 + n*(-1 + gamma_0)))/(2.*brel*(1 + n)*gamma_0*xi*theta) + (brel*pow(E,lambda)*(-1 + n*(-1 + gamma_0))*xi*pow(theta,n))/gamma_0 + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*xi*pow(theta,1 + n);
//----------------------
double alphaV_3 = (pow(E,nu)*(l*alpha_3 + xi*(2 + beta_3)))/((1 + n)*xi*pow(Omega,2));
//----------------------

//----------------------
double alphaH_0 = ((-2*(2 + pow(E,lambda)*(-1 + l))*(2 + l) + ((pow(E,lambda)*(-1 + l) - 2*l)*(2 + l) + l*lambda_der_1*xi)*alpha_0 + lambda_der_1*xi*(8 + xi*beta_0) - 2*xi*((3 + 2*l)*beta_0 + xi*(alphaW_0*beta_1 + alphaV_0*beta_2 + beta_der_0)))/pow(xi,2) - (4*pow(E,lambda)*(1 + n)*pow(theta,-1 + n)*(-1 + brel*(-1 + 2*gamma_0)*theta))/gamma_0)/(2.*beta_3);
//----------------------
double alphaH_1 = (-2*(alphaV_1*beta_2 + beta_der_1) + (((pow(E,lambda)*(-1 + l) - 2*l)*(2 + l) + l*lambda_der_1*xi)*alpha_1 + xi*(-6 - 4*l + lambda_der_1*xi - 2*xi*alphaW_1)*beta_1 + (2*pow(E,lambda/2.)*(1 + n)*pow(theta,-1 + n)*(1 - pow(E,lambda) + brel*(-2*n*theta_der_1*gamma_0*xi + theta - pow(E,lambda)*theta - 2*brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(theta,1 + n)*(1 + brel*theta))))/(brel*gamma_0))/pow(xi,2))/(2.*beta_3);
//----------------------
double alphaH_2 = ((((pow(E,lambda)*(-1 + l) - 2*l)*(2 + l) + l*lambda_der_1*xi)*alpha_2)/pow(xi,2) - 2*alphaW_2*beta_1 + (lambda_der_1 - (2*(3 + 2*l + xi*alphaV_2))/xi)*beta_2 - 2*beta_der_2 - (4*pow(E,lambda - nu)*pow(1 + n,2)*pow(Omega,2)*pow(theta,-1 + n)*(1 + brel*theta))/gamma_0)/(2.*beta_3);
//----------------------
double alphaH_3 = (((pow(E,lambda)*(-1 + l) - 2*l)*(2 + l) + l*lambda_der_1*xi)*alpha_3 + xi*(lambda_der_1*xi*beta_3 - 2*(2 + (3 + 2*l)*beta_3 + xi*(beta_0 + alphaW_3*beta_1 + alphaV_3*beta_2 + beta_der_3))))/(2.*pow(xi,2)*beta_3);
//----------------------

alpha_H2_vals = {alphaH_0,alphaH_1,alphaH_2,alphaH_3};
alpha_V_vals = {alphaV_0, alphaV_1, alphaV_2, alphaV_3};
alpha_W_vals = {alphaW_0, alphaW_1, alphaW_2, alphaW_3 };

}
//======================================================================================================
RHS_PF_SOURCELESS::RHS_PF_SOURCELESS(
    double l, double Omega, 
    double n, double brel,
    double gamma_0,
    std::function <double(double)> thetafunc,  
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc):
l{l},
Omega{Omega},
n{n},
brel{brel},
gamma_0{gamma_0},
thetafunc{thetafunc},
lambdafunc{lambdafunc},
nufunc{nufunc}
{
}
//-------------------------------------------------
RHS_PF_SOURCELESS::~RHS_PF_SOURCELESS(void)
{
    
}
//---------------------------------------------------

//---------------------------------------------------
void RHS_PF_SOURCELESS::get_dydx(double xi, const std::vector<double> &y, std::vector<double> &dydx)
{   
    double theta = thetafunc(xi);
    double lambda = lambdafunc(xi);
    double nu = nufunc(xi);

    vector<double> alpha_H2_vals(4), alpha_W_vals(4), alpha_V_vals(4);

    alpha_H2_W_V_vals_source_less(l, Omega, xi, theta, n, brel, gamma_0, lambda, nu, alpha_H2_vals, alpha_W_vals, alpha_V_vals);

    dydx[0] = y[3];
    dydx[1] = alpha_W_vals[0]*y[0] + alpha_W_vals[1]*y[1] + alpha_W_vals[2]*y[2] + alpha_W_vals[3]*y[3];
    dydx[2] = alpha_V_vals[0]*y[0] + alpha_V_vals[1]*y[1] + alpha_V_vals[2]*y[2] + alpha_V_vals[3]*y[3];
    dydx[3] = alpha_H2_vals[0]*y[0] + alpha_H2_vals[1]*y[1] + alpha_H2_vals[2]*y[2] + alpha_H2_vals[3]*y[3];

}
// //---------------------------------------------------
void RHS_PF_SOURCELESS::operator()(const std::vector<double> &y , std::vector<double> &dydt , double r )
{
    get_dydx(r, y, dydt);
}
//---------------------------------------------------
//==========================================================================
//-----------------Source and functions for bulk viscosity------------------------
void alpha_funcs_source_bulk(double l, double Omega, double xi, double theta, 
 double n, double brel, double lambda, double nu, double &alpha_4, double &alpha_7)
{

double nu_der_1 = (-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*pow(xi,2)*pow(theta,1 + n) + 2*pow(brel,2)*pow(E,lambda)*n*pow(xi,2)*pow(theta,1 + n))/xi;

double denom = 
pow(E,nu)*(pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3)) + 4*brel*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*pow(E,lambda)*(-1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*(1 + n)*nu_der_1*pow(xi,3)*pow(Omega,2) + brel*(1 + n)*pow(nu_der_1,2)*pow(xi,4)*pow(Omega,2) + 4*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*pow(xi,4)*pow(Omega,4))
;
//----------------------
 //----------------------
double alpha_val_4 = 16*pow(E,nu)*Pi*(pow(E,nu)*l*(1 + l) - 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2));
alpha_4=alpha_val_4/denom ;
//----------------------
double alpha_val_7 = -16*pow(E,lambda + nu)*Pi*(pow(E,nu)*l*(1 + l) - 2*brel*(1 + n)*pow(xi,2)*pow(Omega,2));
alpha_7=alpha_val_7/denom ;
//----------------------
}

void beta_funcs_and_ders_source_bulk(double l, double Omega, double xi, double theta, 
 double n, double brel, double lambda, double nu, 
double &beta_4, double &beta_7,
double &beta_der_4, double &beta_der_7)
{
//-----------------------------------------------
double theta_der_1 = -0.5*((1 + brel*theta)*(-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(theta,1 + n)))/(brel*(1 + n)*xi);

double lambda_der_1 = (1 - pow(E,lambda) + 2*brel*pow(E,lambda)*pow(xi,2)*pow(theta,n) + 2*brel*pow(E,lambda)*n*pow(xi,2)*pow(theta,n))/xi;

double nu_der_1 = (-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*pow(xi,2)*pow(theta,1 + n) + 2*pow(brel,2)*pow(E,lambda)*n*pow(xi,2)*pow(theta,1 + n))/xi;


double nu_der_2 = pow(xi,-2) - pow(E,2*lambda)/pow(xi,2) + brel*pow(E,lambda)*(1 + pow(E,lambda))*(1 + n)*pow(theta,n) - pow(brel,2)*pow(E,lambda)*(-5 + 3*pow(E,lambda))*(1 + n)*pow(theta,1 + n) + 2*pow(brel,3)*pow(E,2*lambda)*pow(1 + n,2)*pow(xi,2)*pow(theta,1 + 2*n) - 2*pow(brel,4)*pow(E,2*lambda)*pow(1 + n,2)*pow(xi,2)*pow(theta,2 + 2*n)
; 
//-----------------------------
double denom = pow(E,nu)*xi*(pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3)) + 4*brel*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*pow(E,lambda)*(-1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*(1 + n)*nu_der_1*pow(xi,3)*pow(Omega,2) + brel*(1 + n)*pow(nu_der_1,2)*pow(xi,4)*pow(Omega,2) + 4*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*pow(xi,4)*pow(Omega,4));

double denom_der = pow(E,nu)*(pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3)) + pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3))*nu_der_1*xi + pow(E,lambda + nu)*l*(-2 - l + 2*pow(l,2) + pow(l,3))*(lambda_der_1 + nu_der_1)*xi + 12*brel*(1 + n)*pow(xi,2)*pow(Omega,2) - 12*brel*pow(E,lambda)*(-1 + l + pow(l,2))*(1 + n)*pow(xi,2)*pow(Omega,2) - 4*brel*pow(E,lambda)*(-1 + l + pow(l,2))*lambda_der_1*(1 + n)*pow(xi,3)*pow(Omega,2) - 12*brel*(1 + n)*nu_der_1*pow(xi,3)*pow(Omega,2) - 4*brel*pow(E,lambda)*(-1 + l + pow(l,2))*(1 + n)*nu_der_1*pow(xi,3)*pow(Omega,2) + brel*(1 + n)*pow(nu_der_1,2)*pow(xi,4)*pow(Omega,2) - 4*brel*(1 + n)*nu_der_2*pow(xi,4)*pow(Omega,2) + brel*(1 + n)*pow(nu_der_1,3)*pow(xi,5)*pow(Omega,2) + 2*brel*(1 + n)*nu_der_1*nu_der_2*pow(xi,5)*pow(Omega,2) + 20*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*pow(xi,4)*pow(Omega,4) + 4*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*(lambda_der_1 - nu_der_1)*pow(xi,5)*pow(Omega,4) + 4*pow(brel,2)*pow(E,lambda - nu)*pow(1 + n,2)*nu_der_1*pow(xi,5)*pow(Omega,4))
;
//----------------------
//----------------------
double beta_val_4 = -16*pow(E,nu)*Pi*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n));
beta_4=beta_val_4/denom ;
//----------------------
double beta_val_7 = 16*pow(E,lambda + nu)*Pi*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n));
beta_7=beta_val_7/denom ;
//----------------------   

//----------------------
double beta_val_der_4 = 16*pow(E,nu)*Pi*(-2*pow(E,nu)*pow(l,2)*(1 + l)*nu_der_1 - 2*brel*pow(E,lambda)*(1 + n)*xi*pow(Omega,2) + 2*brel*(3 + 2*l)*(1 + n)*xi*pow(Omega,2) - brel*pow(E,lambda)*lambda_der_1*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*pow(E,lambda)*(1 + n)*nu_der_1*pow(xi,2)*pow(Omega,2) + brel*(3 + 2*l)*(1 + n)*nu_der_1*pow(xi,2)*pow(Omega,2) - 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,3)*theta_der_1*pow(xi,4)*pow(Omega,2)*pow(theta,n) - 8*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,3)*pow(Omega,2)*pow(theta,1 + n) - 2*pow(brel,3)*pow(E,lambda)*lambda_der_1*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n) - 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*nu_der_1*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n));
beta_der_4=beta_val_der_4/denom - (denom_der/(denom*denom))*beta_val_4 ; 
//----------------------
double beta_val_der_7 = 16*pow(E,lambda + nu)*Pi*(pow(E,nu)*pow(l,2)*(1 + l)*nu_der_1 + 2*brel*pow(E,lambda)*(1 + n)*xi*pow(Omega,2) - 2*brel*(3 + 2*l)*(1 + n)*xi*pow(Omega,2) + brel*pow(E,lambda)*lambda_der_1*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,3)*theta_der_1*pow(xi,4)*pow(Omega,2)*pow(theta,n) + 8*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,3)*pow(Omega,2)*pow(theta,1 + n) + 2*pow(brel,3)*pow(E,lambda)*lambda_der_1*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n) + (lambda_der_1 + nu_der_1)*(pow(E,nu)*pow(l,2)*(1 + l) + brel*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(Omega,2) - brel*(3 + 2*l)*(1 + n)*pow(xi,2)*pow(Omega,2) + 2*pow(brel,3)*pow(E,lambda)*pow(1 + n,2)*pow(xi,4)*pow(Omega,2)*pow(theta,1 + n)));
beta_der_7=beta_val_der_7/denom - (denom_der/(denom*denom))*beta_val_7 ; 
//----------------------
}

void alpha_H2_W_V_vals_source_bulk(double l, double Omega, double xi,
double theta, 
double n, double brel, 
double gamma_0,
double lambda, double nu,
double &alphaH_4, double &alphaH_7, double &alphaH_8, 
double &alphaW_4, double &alphaW_7,  
double &alphaV_4, double &alphaV_7, double &alphaV_8 )
{

vector<double>  beta_vals(4), beta_vals_der(4);
beta_funcs_and_ders_source_less(l, Omega,  xi, theta, n, brel,lambda, nu, beta_vals, beta_vals_der);
double  beta_1 = beta_vals[1], beta_2 = beta_vals[2], beta_3 = beta_vals[3];
//------------------------
double alpha_4 = 0., alpha_7 = 0.;
double beta_4 = 0., beta_7 = 0.;
double beta_der_4 = 0., beta_der_7 = 0.;
//------------------------
double lambda_der_1 = (1 - pow(E,lambda) + 2*brel*pow(E,lambda)*pow(xi,2)*pow(theta,n) + 2*brel*pow(E,lambda)*n*pow(xi,2)*pow(theta,n))/xi;

//-------------------------
alpha_funcs_source_bulk(l, Omega, xi, theta, n, brel, lambda, nu, alpha_4, alpha_7);

beta_funcs_and_ders_source_bulk(l, Omega, xi, theta, n, brel, lambda, nu, beta_4, beta_7, beta_der_4, beta_der_7);
//----------------------
alphaW_4 = brel*pow(E,lambda/2.)*xi*alpha_4;
//----------------------
alphaW_7 = (pow(E,lambda/2.)*(pow(brel,2)*pow(xi,2)*alpha_7 + (4*Pi*pow(theta,-1 - n))/((1 + n)*gamma_0*(1 + brel*theta))))/(brel*xi);
//----------------------

//----------------------
alphaV_4 = (pow(E,nu)*((2*(-3 + 3*pow(E,lambda) + 2*l)*Pi)/(brel*pow(E,lambda)*pow(theta,n)) + (1 + n)*pow(xi,2)*(-8*Pi + l*alpha_4 + xi*beta_4 + brel*(4*Pi + l*alpha_4 + xi*beta_4)*theta)))/(pow(1 + n,2)*pow(xi,3)*pow(Omega,2)*(1 + brel*theta));
//----------------------
alphaV_7 = (pow(E,nu)*pow(theta,-1 - n)*(-2*Pi*(-1 + pow(E,lambda) + brel*(-1 + pow(E,lambda) + 2*l)*gamma_0*theta) + pow(brel,2)*(1 + n)*pow(xi,2)*pow(theta,1 + n)*(gamma_0*(l*alpha_7 + xi*beta_7)*(1 + brel*theta) - 4*pow(E,lambda)*Pi*(1 + brel*gamma_0*theta))))/(pow(brel,2)*pow(1 + n,2)*pow(xi,3)*pow(Omega,2)*(gamma_0 + brel*gamma_0*theta));
//----------------------
alphaV_8 = (4*pow(E,-lambda + nu)*Pi)/(brel*pow(1 + n,2)*pow(xi,2)*pow(Omega,2)*pow(theta,n)*(1 + brel*theta));
//----------------------

//----------------------
alphaH_4 = (((pow(E,lambda)*(-1 + l) - 2*l)*(2 + l) + l*lambda_der_1*xi)*alpha_4 + xi*(lambda_der_1*xi*beta_4 - 2*((3 + 2*l)*beta_4 + xi*(alphaW_4*beta_1 + alphaV_4*beta_2 + beta_der_4))))/(2.*pow(xi,2)*beta_3);
//----------------------
alphaH_7 = (brel*gamma_0*(lambda_der_1*xi*(l*alpha_7 + xi*beta_7) - 2*(l*(2 + l)*alpha_7 + xi*((3 + 2*l)*beta_7 + xi*(alphaW_7*beta_1 + alphaV_7*beta_2 + beta_der_7))))*theta + pow(E,lambda)*(-16*Pi + brel*(-2 + l + pow(l,2))*gamma_0*alpha_7*theta))/(2.*brel*gamma_0*pow(xi,2)*beta_3*theta);
//----------------------
alphaH_8 = -((alphaV_8*beta_2 + beta_4)/beta_3);
//----------------------


}

//-------------------------------------------------------------------------------------------
double S_Omega_bulk(double Omega, double xi,
double theta, 
double n, double brel, 
double gamma_0,
double lambda, double nu,
double fzeta,
double H0,
double V,
double W
)
{

return
-0.25*(brel*(1 + n)*pow(xi,2)*Omega*fzeta*H0)/(pow(E,nu/2.)*Pi*gamma_0*theta) + (brel*pow(1 + n,2)*pow(xi,2)*pow(Omega,3)*fzeta*V)/(4.*pow(E,(3*nu)/2.)*Pi*gamma_0*theta) + (pow(E,-0.5*lambda - nu/2.)*(1 + n)*Omega*fzeta*W*(-1 + pow(E,lambda) + 2*pow(brel,2)*pow(E,lambda)*(1 + n)*pow(xi,2)*pow(theta,1 + n)))/(8.*Pi*gamma_0*theta)
;
}

//===============================================================================================================================================
RHS_PF_BULK::RHS_PF_BULK(
    double l, double Omega, 
    double n, double brel,
    double gamma_0,
    std::function <double(double)> thetafunc,  
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc,
    std::function <double(double, double)> fzetafunc,
    std::function <double(double)> H00func,
    std::function <double(double)> W0func,
    std::function <double(double)> V0func):
l{l},
Omega{Omega},
n{n},
brel{brel},
gamma_0{gamma_0},
thetafunc{thetafunc},
lambdafunc{lambdafunc},
nufunc{nufunc},
fzetafunc{fzetafunc},
H00func{H00func},
W0func{W0func},
V0func{V0func}
{
}
//-------------------------------------------------
RHS_PF_BULK::~RHS_PF_BULK(void)
{
    
}
//---------------------------------------------------
//---------------------------------------------------
void RHS_PF_BULK::get_dydx(double xi, const std::vector<double> &y, std::vector<double> &dydx)
{   
    double theta = thetafunc(xi);
    double lambda = lambdafunc(xi);
    double nu = nufunc(xi);
    double lambda_der_1 = (1 - pow(E,lambda) + 2*brel*pow(E,lambda)*pow(xi,2)*pow(theta,n) + 2*brel*pow(E,lambda)*n*pow(xi,2)*pow(theta,n))/xi;


    vector<double> alpha_H2_vals(4), alpha_W_vals(4), alpha_V_vals(4);

    alpha_H2_W_V_vals_source_less(l, Omega, xi, theta, n, brel, gamma_0, lambda, nu, alpha_H2_vals, alpha_W_vals, alpha_V_vals);

    dydx[0] = y[3];
    dydx[1] = alpha_W_vals[0]*y[0] + alpha_W_vals[1]*y[1] + alpha_W_vals[2]*y[2] + alpha_W_vals[3]*y[3];
    dydx[2] = alpha_V_vals[0]*y[0] + alpha_V_vals[1]*y[1] + alpha_V_vals[2]*y[2] + alpha_V_vals[3]*y[3];
    dydx[3] = alpha_H2_vals[0]*y[0] + alpha_H2_vals[1]*y[1] + alpha_H2_vals[2]*y[2] + alpha_H2_vals[3]*y[3];

//------------------------------------------------------------------
    double alphaH_4 =0., alphaH_7 =0., alphaH_8 =0., alphaW_4 =0., alphaW_7 =0.,  alphaV_4 =0., alphaV_7 =0., alphaV_8 =0.;

    alpha_H2_W_V_vals_source_bulk(l, Omega, xi, theta, n, brel, gamma_0, lambda, nu, alphaH_4, alphaH_7, alphaH_8, 
    alphaW_4, alphaW_7, alphaV_4, alphaV_7, alphaV_8 );
//------------------------------------------------------------------
    double H00 = H00func(xi);
    double V0 = V0func(xi);
    double W0 = W0func(xi);

    double fzeta = fzetafunc(theta,n);

    double S_Omega =  S_Omega_bulk(Omega, xi, theta, n, brel, gamma_0, lambda, nu, fzeta, H00, V0, W0);
    double S0 = S_Omega*pow(E,lambda);
//-------------------------------------------------------------------
    dydx[1] += alphaW_4*S0 + alphaW_7*S_Omega;

    dydx[2] += alphaV_4*S0 + alphaV_7*S_Omega + alphaV_8*S_Omega*pow(E,lambda)*lambda_der_1;

    dydx[3] += alphaH_4*S0 + alphaH_7*S_Omega + alphaH_8*S_Omega*pow(E,lambda)*lambda_der_1;



}
// //---------------------------------------------------
void RHS_PF_BULK::operator()(const std::vector<double> &y , std::vector<double> &dydt , double r )
{
    get_dydx(r, y, dydt);
}
//---------------------------------------------------