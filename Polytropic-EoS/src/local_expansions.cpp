#include <iostream>
using std::cout;
using std::endl;
#include <cmath>
using std::pow;
using std::abs;
using std::log;
#include<cassert>
#include <Eigen/Dense>
#include "Li2.hpp"
#include "Li3.hpp"
// #include "Li2.cpp"
// #include "Li3.cpp"
using namespace polylogarithm;
const double E = 2.718281828459045;
const double Pi = 3.14159265358979323846;
//-------------------------------------------------
#include "local_expansions.hpp"
//-------------------------------------------------
//-------------------------------------------------
//-------------------------------------------------
void get_matrix_ord_2_origin(double l, double Omega, double nu0, double n, double brel, 
double gamma_0,
double xi, Eigen::MatrixXd &A)
{
 //------------------
A(0,0) = 1 - ((1 + n)*pow(xi,2)*(pow(E,nu0)*(1 + l)*(3 + 3*pow(brel,2)*(9 + l)*gamma_0 + brel*(3 + (15 - 3*l - 2*pow(l,2))*gamma_0)) + 3*brel*(9 + l)*gamma_0*pow(Omega,2)))/(6.*pow(E,nu0)*(1 + l)*(3 + 2*l)*gamma_0);
//------------------
A(0,1) = ((-1 - brel)*(1 + n)*pow(xi,2)*((1 + 3*brel)*pow(E,nu0)*l*(1 + l)*(-1 + n*(-1 + gamma_0)) - 3*(1 + n)*(-1 + 9*brel*gamma_0 + l*(-1 + brel*gamma_0))*pow(Omega,2)))/(6.*pow(E,nu0)*l*(1 + l)*(3 + 2*l)*gamma_0);
//------------------
A(1,0) = (pow(xi,2)*(3*(-2 - 6*brel*gamma_0 + l*(-1 + brel*gamma_0)) + ((1 + 3*brel)*pow(E,nu0)*l*(1 + l)*(-1 + n*(-1 + gamma_0)))/((1 + n)*pow(Omega,2))))/(6.*(3 + 2*l)*gamma_0);
//------------------
A(1,1) = (pow(1 + 3*brel,2)*pow(E,nu0)*pow(l,2)*(1 + l)*(1 + n - n*gamma_0)*pow(xi,2) + 18*l*(3 + 2*l)*gamma_0*pow(Omega,2) + 3*l*(1 + n + n*gamma_0 + l*n*gamma_0 + 6*pow(brel,2)*(-1 + l)*(1 + n)*gamma_0 + brel*(3 + 2*(-4 + 3*l + pow(l,2))*gamma_0 + n*(3 + (-5 + 9*l + 2*pow(l,2))*gamma_0)))*pow(xi,2)*pow(Omega,2) - (9*(2 + l)*(1 + n)*pow(xi,2)*pow(Omega,4))/pow(E,nu0))/(18.*l*(3 + 2*l)*gamma_0*pow(Omega,2));
//------------------
A(2,0) = (pow(xi,2)*(3*(1 + l - 9*brel*gamma_0 - brel*l*gamma_0) + ((-1 - 3*brel)*pow(E,nu0)*(3 + 4*l + pow(l,2))*(-1 + n*(-1 + gamma_0)))/((1 + n)*pow(Omega,2))))/(6.*(1 + l)*(3 + 2*l)*gamma_0);
//------------------
A(2,1) = (pow(E,nu0)*l*(3 + 4*l + pow(l,2))*(-1 + n*(-1 + gamma_0))*pow(xi + 3*brel*xi,2) - 18*(1 + l)*(3 + 2*l)*gamma_0*pow(Omega,2) - 3*(6*pow(brel,2)*(-3 + 2*l + pow(l,2))*(1 + n)*gamma_0 + (1 + l)*(-3 + n*(-3 + (3 + l)*gamma_0)) + brel*(2*pow(l,3)*(1 + n)*gamma_0 + pow(l,2)*(8 + 11*n)*gamma_0 + 3*l*(-3 - 3*n + 4*n*gamma_0) - 3*(3 + 3*n + 10*gamma_0 + 7*n*gamma_0)))*pow(xi,2)*pow(Omega,2) + (9*(1 + l)*(1 + n)*pow(xi,2)*pow(Omega,4))/pow(E,nu0))/(18.*l*(1 + l)*(3 + 2*l)*gamma_0*pow(Omega,2));
//------------------
A(3,0) = -0.3333333333333333*((1 + n)*xi*(pow(E,nu0)*(1 + l)*(3 + 3*pow(brel,2)*(9 + l)*gamma_0 + brel*(3 + (15 - 3*l - 2*pow(l,2))*gamma_0)) + 3*brel*(9 + l)*gamma_0*pow(Omega,2)))/(pow(E,nu0)*(1 + l)*(3 + 2*l)*gamma_0);
//------------------
A(3,1) = ((-1 - brel)*(1 + n)*xi*((1 + 3*brel)*pow(E,nu0)*l*(1 + l)*(-1 + n*(-1 + gamma_0)) - 3*(1 + n)*(-1 + 9*brel*gamma_0 + l*(-1 + brel*gamma_0))*pow(Omega,2)))/(3.*pow(E,nu0)*l*(1 + l)*(3 + 2*l)*gamma_0);
//------------------
}

void get_matrix_ord_2_origin_just_H_W(double l, double Omega, double nu0, double n, double brel, 
double gamma_0,
double xi, Eigen::MatrixXd &A)
{
//------------------
A(0,0) = 1 - ((1 + n)*pow(xi,2)*(pow(E,nu0)*(1 + l)*(3 + 3*pow(brel,2)*(9 + l)*gamma_0 + brel*(3 + (15 - 3*l - 2*pow(l,2))*gamma_0)) + 3*brel*(9 + l)*gamma_0*pow(Omega,2)))/(6.*pow(E,nu0)*(1 + l)*(3 + 2*l)*gamma_0);
//------------------
A(0,1) = ((-1 - brel)*(1 + n)*pow(xi,2)*((1 + 3*brel)*pow(E,nu0)*l*(1 + l)*(-1 + n*(-1 + gamma_0)) - 3*(1 + n)*(-1 + 9*brel*gamma_0 + l*(-1 + brel*gamma_0))*pow(Omega,2)))/(6.*pow(E,nu0)*l*(1 + l)*(3 + 2*l)*gamma_0);
//------------------
A(1,0) = (pow(xi,2)*(3*(-2 - 6*brel*gamma_0 + l*(-1 + brel*gamma_0)) + ((1 + 3*brel)*pow(E,nu0)*l*(1 + l)*(-1 + n*(-1 + gamma_0)))/((1 + n)*pow(Omega,2))))/(6.*(3 + 2*l)*gamma_0);
//------------------
A(1,1) = (pow(1 + 3*brel,2)*pow(E,nu0)*pow(l,2)*(1 + l)*(1 + n - n*gamma_0)*pow(xi,2) + 18*l*(3 + 2*l)*gamma_0*pow(Omega,2) + 3*l*(1 + n + n*gamma_0 + l*n*gamma_0 + 6*pow(brel,2)*(-1 + l)*(1 + n)*gamma_0 + brel*(3 + 2*(-4 + 3*l + pow(l,2))*gamma_0 + n*(3 + (-5 + 9*l + 2*pow(l,2))*gamma_0)))*pow(xi,2)*pow(Omega,2) - (9*(2 + l)*(1 + n)*pow(xi,2)*pow(Omega,4))/pow(E,nu0))/(18.*l*(3 + 2*l)*gamma_0*pow(Omega,2));
//------------------
}
//===================================================================================================================
SERIES_ORIGIN_SOURCELESS::SERIES_ORIGIN_SOURCELESS(double l, double Omega, double nu0, double n, double brel, double gamma_0):
l{l},
Omega{Omega},
nu0{nu0},
n{n},
brel{brel},
gamma_0{gamma_0}
{

}

//-------------------------------------------------
SERIES_ORIGIN_SOURCELESS::~SERIES_ORIGIN_SOURCELESS(void)
{
}
//-------------------------------------------------
void SERIES_ORIGIN_SOURCELESS::operator()( double xi, const std::vector<double> yin, std::vector<double> &y )
{   
    Eigen::MatrixXd A(4,2);

    get_matrix_ord_2_origin(l, Omega, nu0, n, brel, gamma_0, xi, A);

    for(int i=0; i<4; ++i)
    {
        y[i] = A(i,0)*yin[0] + A(i,1)*yin[1];
    }
    
}
//-------------------------------------------------
SERIES_SURF::SERIES_SURF(double l, double Omega, double n, double brel, double gamma_0, double muR0, double xi1):
l{l},
Omega{Omega},
n{n},
brel{brel},
gamma_0{gamma_0},
muR0{muR0},
xi1{xi1}
{

}

//-------------------------------------------------
SERIES_SURF::~SERIES_SURF(void)
{
}
//-------------------------------------------------
void SERIES_SURF::operator()( double xi, const std::vector<double> yin, std::vector<double> &y )
{   

//-----------------------------
double H_R0= yin[0];
double H_R1= yin[1];
double W_R0= yin[2];
//-----------------------------

//-------------------------
double x1 = xi1- 2*brel*muR0*(1 + n);
double thetaR1 = -(muR0/((2*brel*muR0*(1 + n) - xi1)*xi1));
// double thetaR2 = -0.5*(muR0*(brel*muR0 + 2*brel*muR0*n - 2*xi1))/(pow(2*brel*muR0 + 2*brel*muR0*n - xi1,2)*pow(xi1,2));
//-----------------------------
double a0 = (4*pow(muR0,2)*Pi*sqrt((1 + n)*x1)*pow(2*brel*muR0*(1 + n) + x1,-5 + l)*(2*brel*muR0*pow(1 + n,2.5)*sqrt(x1) - (1 + n)*sqrt((1 + n)*x1)*(2*brel*muR0*(1 + n) + x1))*(1 + n - n*gamma_0))/(brel*pow(1 + n,2)*thetaR1*pow(x1,2)*sqrt(x1/(2*brel*muR0*(1 + n) + x1))*gamma_0) ;

//------------------------------------------------------
double b0 = (pow(x1/(2*brel*muR0*(1 + n) + x1),1.5)*(1 + n + gamma_0))/(x1*gamma_0);

//------------------------------------------------------
double a1 = (4*Pi*sqrt((1 + n)*x1)*pow(2*brel*muR0*(1 + n) + x1,-1 + l)*(thetaR1*sqrt((1 + n)*x1)*pow(2*brel*muR0*(1 + n) + x1,2)*gamma_0 + muR0*(-sqrt((1 + n)*x1) - n*sqrt((1 + n)*x1) + n*sqrt((1 + n)*x1)*gamma_0 - 2*brel*sqrt(1 + n)*thetaR1*sqrt(x1)*(2*brel*muR0*(1 + n) + x1)*gamma_0 - 2*brel*n*sqrt(1 + n)*thetaR1*sqrt(x1)*(2*brel*muR0*(1 + n) + x1)*gamma_0))*pow(Omega,2))/(brel*(1 + n)*thetaR1*pow(x1,2)*gamma_0) ;

//------------------------------------------------------
double b1 = ((1 + n)*pow(2*brel*muR0*(1 + n) + x1,2)*pow(Omega,2))/(muR0*gamma_0);

//------------------------------------------------------
double s0 = (-4*Pi*W_R0*sqrt((1 + n)*x1)*sqrt(x1/(2*brel*muR0*(1 + n) + x1))*pow(2*brel*muR0*(1 + n) + x1,-5 + l)*(7*brel*pow(muR0,3)*pow(1 + n,2.5)*sqrt(x1)*(-1 + n*(-1 + gamma_0)) - 3*pow(muR0,2)*(1 + n)*sqrt((1 + n)*x1)*(2*brel*muR0*(1 + n) + x1)*(-1 + n*(-1 + gamma_0)) - 4*pow(brel,2)*(1 + l)*pow(muR0,3)*pow(1 + n,3.5)*thetaR1*sqrt(x1)*(2*brel*muR0*(1 + n) + x1)*gamma_0 + 2*brel*(1 + 2*l)*pow(muR0,2)*pow(1 + n,2.5)*thetaR1*sqrt(x1)*pow(2*brel*muR0*(1 + n) + x1,2)*gamma_0 - l*muR0*(1 + n)*thetaR1*sqrt((1 + n)*x1)*pow(2*brel*muR0*(1 + n) + x1,3)*gamma_0 - 2*brel*muR0*pow(1 + n,2.5)*thetaR1*sqrt(x1)*pow(2*brel*muR0*(1 + n) + x1,5)*gamma_0*pow(Omega,2) + (1 + n)*thetaR1*sqrt((1 + n)*x1)*pow(2*brel*muR0*(1 + n) + x1,6)*gamma_0*pow(Omega,2)))/(brel*pow(1 + n,2)*thetaR1*pow(x1,3)*gamma_0) + (4*H_R1*Pi*sqrt((1 + n)*x1)*pow(2*brel*muR0*(1 + n) + x1,-2 + l)*((-2*pow(brel,1.5)*pow(x1,2.5)*(2*brel*muR0*(1 + n) + x1))/pow(brel*(1 + n),1.5) - (muR0*(2*brel*muR0*pow(1 + n,1.5)*sqrt(x1) - sqrt((1 + n)*x1)*(2*brel*muR0*(1 + n) + x1))*(1 + n - n*gamma_0))/(pow(1 + n,2)*thetaR1*gamma_0) + (4*pow(brel,1.5)*l*pow(x1,4.5)*(l*muR0 + pow(l,2)*muR0 - 2*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)))/(sqrt(brel*(1 + n))*(l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(2*brel*muR0*(1 + n) + x1,2) + 36*pow(brel,3)*pow(muR0,2)*pow(1 + n,3)*pow(2*brel*muR0*(1 + n) + x1,2)*pow(Omega,2) - 4*brel*(-2 + l + pow(l,2))*(1 + n)*(2*brel*muR0*(1 + n) + x1)*(l*muR0 + pow(l,2)*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + 4*pow(brel,2)*pow(1 + n,2)*(2*pow(l,3)*pow(muR0,2) + pow(l,4)*pow(muR0,2) - 8*muR0*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2) + pow(2*brel*muR0*(1 + n) + x1,6)*pow(Omega,4) - pow(l,2)*muR0*(muR0 - 2*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) - 2*l*muR0*(muR0 - pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2))))) + (2*pow(brel,1.5)*pow(x1,3.5)*(l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(2*brel*muR0*(1 + n) + x1,2) + 4*pow(brel,2)*muR0*pow(1 + n,2)*(pow(l,2)*muR0 + pow(l,3)*muR0 - 3*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2) - 2*l*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) - 2*brel*(1 + l)*(1 + n)*(2*brel*muR0*(1 + n) + x1)*(-2*l*muR0 + 2*pow(l,2)*muR0 + pow(l,3)*muR0 - 2*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2) + l*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2))))/(pow(brel*(1 + n),1.5)*(l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(2*brel*muR0*(1 + n) + x1,2) + 36*pow(brel,3)*pow(muR0,2)*pow(1 + n,3)*pow(2*brel*muR0*(1 + n) + x1,2)*pow(Omega,2) - 4*brel*(-2 + l + pow(l,2))*(1 + n)*(2*brel*muR0*(1 + n) + x1)*(l*muR0 + pow(l,2)*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + 4*pow(brel,2)*pow(1 + n,2)*(2*pow(l,3)*pow(muR0,2) + pow(l,4)*pow(muR0,2) - 8*muR0*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2) + pow(2*brel*muR0*(1 + n) + x1,6)*pow(Omega,4) - pow(l,2)*muR0*(muR0 - 2*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) - 2*l*muR0*(muR0 - pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)))))))/(brel*pow(x1,2)) + (4*H_R0*Pi*sqrt((1 + n)*x1)*pow(2*brel*muR0*(1 + n) + x1,l)*(-((l*pow((1 + n)*x1,1.5))/pow(1 + n,3)) - (4*pow(brel,2)*(1 + l)*pow(muR0,2)*sqrt(1 + n)*pow(x1,1.5))/pow(2*brel*muR0*(1 + n) + x1,2) + (2*brel*(1 + 2*l)*muR0*pow(x1,1.5))/(sqrt(1 + n)*(2*brel*muR0*(1 + n) + x1)) + (2*pow(brel,1.5)*pow(x1,2.5)*(2*brel*muR0*(1 + n) + l*x1))/(pow(brel*(1 + n),1.5)*pow(2*brel*muR0*(1 + n) + x1,2)) - (2*brel*pow(muR0,2)*pow(x1,1.5)*(1 + n - n*gamma_0))/(sqrt(1 + n)*thetaR1*pow(2*brel*muR0*(1 + n) + x1,3)*gamma_0) - (2*pow(brel,1.5)*l*pow(x1,4.5)*((-2*brel*(2 + l + pow(l,2))*muR0*(1 + n) - (-2 + 3*l + pow(l,2))*x1)*(-(l*(1 + l)*x1) + 2*brel*(1 + n)*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + 2*(brel*muR0*(1 + n) - x1)*(pow(l,2)*(1 + l)*(2*brel*muR0*(1 + n) + x1) - 2*brel*(1 + n)*(-(l*muR0) + pow(l,3)*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)))))/(pow(brel*(1 + n),1.5)*pow(2*brel*muR0*(1 + n) + x1,2)*(2*(brel*muR0*(1 + n) - x1)*(pow(l,2)*(1 + l)*pow(x1,2) + 2*brel*(1 + n)*pow(2*brel*muR0*(1 + n) + x1,3)*(brel*muR0*(1 + n) - (1 + l)*x1)*pow(Omega,2)) + (l*(1 + l)*x1 - 2*brel*(1 + n)*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2))*(12*pow(brel,2)*l*pow(muR0,2)*pow(1 + n,2) + (-2 + 3*l + pow(l,2))*pow(2*brel*muR0*(1 + n) + x1,2) - 2*brel*(1 + n)*(2*brel*muR0*(1 + n) + x1)*((-2 + 6*l + pow(l,2))*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2))))) + (4*pow(brel,1.5)*pow(x1,3.5)*(l*pow(2*brel*muR0*(1 + n) + x1,2)*((2 + 3*l + pow(l,2))*muR0 - 2*(-1 + l)*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + 4*pow(brel,2)*(1 + l)*pow(muR0,2)*pow(1 + n,2)*(pow(l,2)*muR0 + pow(l,3)*muR0 + 6*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2) - 2*l*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) - brel*(1 + n)*(2*brel*muR0*(1 + n) + x1)*(6*pow(l,3)*pow(muR0,2) + 2*pow(l,4)*pow(muR0,2) + pow(l,2)*muR0*(8*muR0 - 11*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + 2*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)*(5*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + l*muR0*(4*muR0 + 9*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)))))/(sqrt(brel*(1 + n))*pow(2*brel*muR0*(1 + n) + x1,3)*(l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(2*brel*muR0*(1 + n) + x1,2) + 36*pow(brel,3)*pow(muR0,2)*pow(1 + n,3)*pow(2*brel*muR0*(1 + n) + x1,2)*pow(Omega,2) - 4*brel*(-2 + l + pow(l,2))*(1 + n)*(2*brel*muR0*(1 + n) + x1)*(l*muR0 + pow(l,2)*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + 4*pow(brel,2)*pow(1 + n,2)*(2*pow(l,3)*pow(muR0,2) + pow(l,4)*pow(muR0,2) - 8*muR0*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2) + pow(2*brel*muR0*(1 + n) + x1,6)*pow(Omega,4) - pow(l,2)*muR0*(muR0 - 2*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) - 2*l*muR0*(muR0 - pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)))))))/(brel*pow(x1,3));

//------------------------------------------------------
double s1 = -((W_R0*(-(pow(l,2)*muR0*x1*gamma_0) + pow(2*brel*muR0*(1 + n) + x1,3)*(brel*muR0*pow(1 + n,2) + x1*(-3 - 3*n + gamma_0))*pow(Omega,2) + l*x1*gamma_0*(8*pow(brel,3)*pow(muR0,3)*pow(1 + n,3)*pow(Omega,2) + 12*pow(brel,2)*pow(muR0,2)*pow(1 + n,2)*x1*pow(Omega,2) + pow(x1,3)*pow(Omega,2) + muR0*(-1 + 6*brel*(1 + n)*pow(x1,2)*pow(Omega,2)))))/(sqrt(x1/(2*brel*muR0*(1 + n) + x1))*pow(2*brel*muR0*(1 + n) + x1,6)*gamma_0*pow(Omega,2))) + H_R1*x1*(-((2*brel*muR0*(1 + n) + x1)/(muR0*gamma_0)) + (4*pow(brel,2)*(1 + n)*x1*(l*muR0 + pow(l,2)*muR0 - 2*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)))/(l*(-2 - l + 2*pow(l,2) + pow(l,3))*pow(2*brel*muR0*(1 + n) + x1,2) + 36*pow(brel,3)*pow(muR0,2)*pow(1 + n,3)*pow(2*brel*muR0*(1 + n) + x1,2)*pow(Omega,2) - 4*brel*(-2 + l + pow(l,2))*(1 + n)*(2*brel*muR0*(1 + n) + x1)*(l*muR0 + pow(l,2)*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + 4*pow(brel,2)*pow(1 + n,2)*(2*pow(l,3)*pow(muR0,2) + pow(l,4)*pow(muR0,2) - 8*muR0*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2) + pow(2*brel*muR0*(1 + n) + x1,6)*pow(Omega,4) - pow(l,2)*muR0*(muR0 - 2*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) - 2*l*muR0*(muR0 - pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2))))) + (H_R0*(2*pow(brel,2)*muR0*(1 + n)*pow(2*brel*muR0*(1 + n) + x1,3) - brel*pow(2*brel*muR0*(1 + n) + x1,4) + (2*brel*(1 + n)*x1*pow(2*brel*muR0*(1 + n) + x1,3))/gamma_0 + (4*pow(brel,2)*l*pow(muR0,2)*(1 + n))/pow(Omega,2) + (4*pow(brel,2)*pow(l,2)*pow(muR0,2)*(1 + n))/pow(Omega,2) + (4*brel*l*muR0*x1)/pow(Omega,2) + (4*brel*pow(l,2)*muR0*x1)/pow(Omega,2) - (l*pow(2*brel*muR0*(1 + n) + x1,2))/((1 + n)*pow(Omega,2)) - (pow(l,2)*pow(2*brel*muR0*(1 + n) + x1,2))/((1 + n)*pow(Omega,2)) - (2*brel*pow(x1,2)*pow(2*brel*muR0*(1 + n) + x1,3)*((-2*brel*(2 + l + pow(l,2))*muR0*(1 + n) - (-2 + 3*l + pow(l,2))*x1)*(-(l*(1 + l)*x1) + 2*brel*(1 + n)*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)) + 2*(brel*muR0*(1 + n) - x1)*(pow(l,2)*(1 + l)*(2*brel*muR0*(1 + n) + x1) - 2*brel*(1 + n)*(-(l*muR0) + pow(l,3)*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2)))))/(2*(brel*muR0*(1 + n) - x1)*(pow(l,2)*(1 + l)*pow(x1,2) + 2*brel*(1 + n)*pow(2*brel*muR0*(1 + n) + x1,3)*(brel*muR0*(1 + n) - (1 + l)*x1)*pow(Omega,2)) + (l*(1 + l)*x1 - 2*brel*(1 + n)*pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2))*(12*pow(brel,2)*l*pow(muR0,2)*pow(1 + n,2) + (-2 + 3*l + pow(l,2))*pow(2*brel*muR0*(1 + n) + x1,2) - 2*brel*(1 + n)*(2*brel*muR0*(1 + n) + x1)*((-2 + 6*l + pow(l,2))*muR0 + pow(2*brel*muR0*(1 + n) + x1,3)*pow(Omega,2))))))/(x1*pow(2*brel*muR0*(1 + n) + x1,3));

//-----------------------------
y[0] = H_R0 + H_R1*(-xi + xi1);
//----------------------
y[1] = W_R0 - ((-(b1*s0) + b0*s1)*(-xi + xi1))/(a1*b0 - a0*b1);
//----------------------
y[2] = -(((a1*s0 - a0*s1)*(-xi + xi1))/(a1*b0 - a0*b1)) + (H_R0*(1/(1 + n) - (2*brel*muR0)/xi1))/pow(Omega,2) - (muR0*W_R0*sqrt(1 - (2*brel*muR0*(1 + n))/xi1))/(pow(xi1,3)*pow(Omega,2));
//----------------------
y[3] = -H_R1;
//----------------------

    
}
//-------------------------------------------------
void obtain_H0_and_H1_ext(double l, double Omega, double n, double brel, double mu1, double xi1, std::vector<double> &H0_ext, std::vector<double> &H1_ext )
{
    assert(abs(l-2)<1e-6); //Exterior problem only implemented for l=2.

    double z = 1 - (2*brel*mu1*(1 + n))/xi1;
    double dzdxi = (2*brel*mu1*(1 + n))/pow(xi1,2);
    double epsilon = mu1*Omega*pow(brel*(1+n),3./2.);
    // double epsilon = 0.0;
 
    double C1 = brel*(n+1)*mu1/xi1;

    double zeta_3 = 1.202056903159594285399738161511449990764;

    double Li_2 = Li2(z);
    double Li_3 = Li3(z);


   H0_ext[0] = 
    z + (pow(epsilon,2)*(8693 - 174290*z + 308389*pow(z,2) + 13440*pow(Pi,2)*pow(z,2) - 208448*pow(z,3) - 26880*pow(Pi,2)*pow(z,3) + 10555*pow(z,4) + 13440*pow(Pi,2)*pow(z,4) + 37010*pow(z,5) - 3029*pow(z,6) - 80640*pow(-1 + z,2)*pow(z,2)*Li_2 + 6720*log(2) - 67200*z*log(2) + 73152*pow(z,2)*log(2) + 82176*pow(z,3)*log(2) - 155328*pow(z,4)*log(2) + 67200*pow(z,5)*log(2) - 6720*pow(z,6)*log(2) + 192*pow(-1 + z,2)*(-35 + 280*z + 214*pow(z,2) - 280*pow(z,3) + 35*pow(z,4))*log(1 - z) + 6720*log(z) - 67200*z*log(z) + 36804*pow(z,2)*log(z) + 154872*pow(z,3)*log(z) - 191676*pow(z,4)*log(z) + 67200*pow(z,5)*log(z) - 6720*pow(z,6)*log(z) - 80640*pow(z,2)*log(2)*log(z) + 161280*pow(z,3)*log(2)*log(z) - 80640*pow(z,4)*log(2)*log(z) - 40320*pow(z,2)*pow(log(z),2) + 80640*pow(z,3)*pow(log(z),2) - 40320*pow(z,4)*pow(log(z),2)))/(20160.*pow(-1 + z,2)*z)
    ;
    
//------------------
    H1_ext[0] = 
    1 + (pow(epsilon,2)*(1973 + 41121*z + 3387*pow(z,2) - 13440*pow(Pi,2)*pow(z,2) - 82713*pow(z,3) + 40320*pow(Pi,2)*pow(z,3) + 232707*pow(z,4) - 40320*pow(Pi,2)*pow(z,4) - 241033*pow(z,5) + 13440*pow(Pi,2)*pow(z,5) + 95885*pow(z,6) - 9087*pow(z,7) - 80640*pow(-1 + z,3)*pow(z,2)*Li_2 + 6720*log(2) - 20160*z*log(2) + 141888*pow(z,2)*log(2) - 479424*pow(z,3)*log(2) + 707904*pow(z,4)*log(2) - 504768*pow(z,5)*log(2) + 168000*pow(z,6)*log(2) - 20160*pow(z,7)*log(2) + 192*pow(-1 + z,3)*(35 + 634*pow(z,2) - 560*pow(z,3) + 105*pow(z,4))*log(1 - z) + 6720*log(z) - 20160*z*log(z) + 178236*pow(z,2)*log(z) - 588468*pow(z,3)*log(z) + 816948*pow(z,4)*log(z) - 541116*pow(z,5)*log(z) + 168000*pow(z,6)*log(z) - 20160*pow(z,7)*log(z) + 80640*pow(z,2)*log(2)*log(z) - 241920*pow(z,3)*log(2)*log(z) + 241920*pow(z,4)*log(2)*log(z) - 80640*pow(z,5)*log(2)*log(z) + 40320*pow(z,2)*pow(log(z),2) - 120960*pow(z,3)*pow(log(z),2) + 120960*pow(z,4)*pow(log(z),2) - 40320*pow(z,5)*pow(log(z),2)))/(20160.*pow(-1 + z,3)*pow(z,2))
    ;


    H1_ext[0] *= dzdxi;

//------------------

    H0_ext[1] = 
    (-5*(-1 + 8*z - 8*pow(z,3) + pow(z,4) + 12*pow(z,2)*log(z)))/(32.*pow(C1,5)*z) + (pow(epsilon,2)*(22050*(-1 + 8*z - 8*pow(z,3) + pow(z,4))*pow(log(z),2) + 88200*pow(z,2)*pow(log(z),3) - 22470*log(2)*(-1 + 8*z - 8*pow(z,3) + pow(z,4) + 12*pow(z,2)*log(z)) + 210*log(1 - z)*(107*(-1 + 8*z - 8*pow(z,3) + pow(z,4)) + 9924*pow(z,2)*log(z) + 2520*pow(z,2)*pow(log(z),2)) + (6*log(z)*(-2450 + 120050*z + 28841*pow(z,2) - 264322*pow(z,3) + 114626*pow(z,4) + 30100*pow(z,5) - 3745*pow(z,6) + 176400*pow(-1 + z,2)*pow(z,2)*Li_2 + 14700*pow(-1 + z,2)*pow(z,2)*(pow(Pi,2) - 6*Li_2 - 6*log(1 - z)*log(z))))/pow(-1 + z,2) + (-31517 - 1060977*z + 1157944*pow(z,2) + 302400*pow(Pi,2)*pow(z,2) + 251444*pow(z,3) - 302400*pow(Pi,2)*pow(z,3) - 509727*pow(z,4) + 54233*pow(z,5) - 1058400*pow(z,2)*zeta_3 + 1058400*pow(z,3)*zeta_3 + 1814400*(-1 + z)*pow(z,2)*Li_2 + 1058400*pow(z,2)*Li_3 - 1058400*pow(z,3)*Li_3 + 210*(35 - 315*z - 1588*pow(z,2) + 2148*pow(z,3) - 315*pow(z,4) + 35*pow(z,5))*(pow(Pi,2) - 6*Li_2 - 6*log(1 - z)*log(z)))/(-1 + z)))/(70560.*pow(C1,5)*z)
    ;

    
 //------------------

    H1_ext[1] = 
    (-5*(1 + 12*pow(z,2) - 16*pow(z,3) + 3*pow(z,4) + 12*pow(z,2)*log(z)))/(32.*pow(C1,5)*pow(z,2)) + (pow(epsilon,2)*(46217 - 7350*pow(Pi,2) - 852021*z + 22050*pow(Pi,2)*z + 931955*pow(z,2) - 200130*pow(Pi,2)*pow(z,2) + 1782909*pow(z,3) + 659190*pow(Pi,2)*pow(z,3) - 4557201*pow(z,4) - 909090*pow(Pi,2)*pow(z,4) + 3906157*pow(z,5) + 597030*pow(Pi,2)*pow(z,5) - 1420715*pow(z,6) - 183750*pow(Pi,2)*pow(z,6) + 162699*pow(z,7) + 22050*pow(Pi,2)*pow(z,7) - 1058400*pow(z,2)*zeta_3 + 3175200*pow(z,3)*zeta_3 - 3175200*pow(z,4)*zeta_3 + 1058400*pow(z,5)*zeta_3 - 1058400*pow(-1 + z,3)*pow(z,2)*Li_3 + 22470*log(2) - 67410*z*log(2) + 337050*pow(z,2)*log(2) - 1190910*pow(z,3)*log(2) + 1954890*pow(z,4)*log(2) - 1550430*pow(z,5)*log(2) + 561750*pow(z,6)*log(2) - 67410*pow(z,7)*log(2) - 22470*log(1 - z) + 67410*z*log(1 - z) - 337050*pow(z,2)*log(1 - z) + 1190910*pow(z,3)*log(1 - z) - 1954890*pow(z,4)*log(1 - z) + 1550430*pow(z,5)*log(1 - z) - 561750*pow(z,6)*log(1 - z) + 67410*pow(z,7)*log(1 - z) + 29400*log(z) - 396900*z*log(z) - 863946*pow(z,2)*log(z) - 88200*pow(Pi,2)*pow(z,2)*log(z) + 2729178*pow(z,3)*log(z) + 264600*pow(Pi,2)*pow(z,3)*log(z) - 2273688*pow(z,4)*log(z) - 264600*pow(Pi,2)*pow(z,4)*log(z) + 136716*pow(z,5)*log(z) + 88200*pow(Pi,2)*pow(z,5)*log(z) + 429450*pow(z,6)*log(z) - 67410*pow(z,7)*log(z) + 269640*pow(z,2)*log(2)*log(z) - 808920*pow(z,3)*log(2)*log(z) + 808920*pow(z,4)*log(2)*log(z) - 269640*pow(z,5)*log(2)*log(z) + 44100*log(1 - z)*log(z) - 132300*z*log(1 - z)*log(z) + 931140*pow(z,2)*log(1 - z)*log(z) - 3146220*pow(z,3)*log(1 - z)*log(z) + 4645620*pow(z,4)*log(1 - z)*log(z) - 3312540*pow(z,5)*log(1 - z)*log(z) + 1102500*pow(z,6)*log(1 - z)*log(z) - 132300*pow(z,7)*log(1 - z)*log(z) - 22050*pow(log(z),2) + 66150*z*pow(log(z),2) - 330750*pow(z,2)*pow(log(z),2) + 1168650*pow(z,3)*pow(log(z),2) - 1918350*pow(z,4)*pow(log(z),2) + 1521450*pow(z,5)*pow(log(z),2) - 551250*pow(z,6)*pow(log(z),2) + 66150*pow(z,7)*pow(log(z),2) - 88200*pow(z,2)*pow(log(z),3) + 264600*pow(z,3)*pow(log(z),3) - 264600*pow(z,4)*pow(log(z),3) + 88200*pow(z,5)*pow(log(z),3) - 1260*pow(-1 + z,3)*Li_2*(35 + 848*pow(z,2) - 560*pow(z,3) + 105*pow(z,4) - 420*pow(z,2)*log(z))))/(70560.*pow(C1,5)*pow(-1 + z,3)*pow(z,2))
    ;
    
    H1_ext[1] *= dzdxi;

/*----------------------------------------*/
//Eric Poisson's potentials
/*----------------------------------------*/
//     H0_ext[0] = 
//     z + (pow(epsilon,2)*(70 - 3430*z + 3989*pow(z,2) - 2074*pow(z,3) + 1538*pow(z,4) - 860*pow(z,5) + 107*pow(z,6) - 2520*pow(-1 + z,2)*pow(z,2)*Li_2 + 6*pow(-1 + z,2)*(-35 + 280*z + 214*pow(z,2) - 280*pow(z,3) + 35*pow(z,4))*log(1 - z) + 210*log(z) - 2100*z*log(z) + 3570*pow(z,2)*log(z) - 3570*pow(z,4)*log(z) + 2100*pow(z,5)*log(z) - 210*pow(z,6)*log(z) - 1260*pow(z,2)*pow(log(z),2) + 2520*pow(z,3)*pow(log(z),2) - 1260*pow(z,4)*pow(log(z),2)))/(630.*pow(-1 + z,2)*z)
//     ;
    
// //------------------
//     H1_ext[0] = 
//     1 + (pow(epsilon,2)*(-140 + 1890*z - 699*pow(z,2) + 1443*pow(z,3) - 3612*pow(z,4) + 4162*pow(z,5) - 2045*pow(z,6) + 321*pow(z,7) - 2520*pow(-1 + z,3)*pow(z,2)*Li_2 + 6*pow(-1 + z,3)*(35 + 634*pow(z,2) - 560*pow(z,3) + 105*pow(z,4))*log(1 - z) + 210*log(z) - 630*z*log(z) + 3150*pow(z,2)*log(z) - 11130*pow(z,3)*log(z) + 18270*pow(z,4)*log(z) - 14490*pow(z,5)*log(z) + 5250*pow(z,6)*log(z) - 630*pow(z,7)*log(z) + 1260*pow(z,2)*pow(log(z),2) - 3780*pow(z,3)*pow(log(z),2) + 3780*pow(z,4)*pow(log(z),2) - 1260*pow(z,5)*pow(log(z),2)))/(630.*pow(-1 + z,3)*pow(z,2))
//     ;


//     H1_ext[0] *= dzdxi;

// //------------------

//     H0_ext[1] = 
//     (-5*(-1 + 8*z - 8*pow(z,3) + pow(z,4) + 12*pow(z,2)*log(z)))/(32.*pow(C1,5)*z) - (pow(epsilon,2)*z*(-6900 - 2568*pow(Pi,2) + (3307 + 8655*z - 20732*pow(z,2) + 18968*pow(z,3) - 7095*pow(z,4) + 857*pow(z,5))/((-1 + z)*pow(z,2)) - 30240*zeta_3 + 30240*Li_3 - (642*(-1 + z)*(1 + z)*(1 - 8*z + pow(z,2))*log(1 - z))/pow(z,2) + (18*(-35 + 280*z + 428*pow(z,2) - 280*pow(z,3) + 35*pow(z,4))*(2*Li_2 + pow(log(1/(1 - z)),2) - pow(log(1 - z),2)))/pow(z,2) + (6*(70 - 3430*z + 3989*pow(z,2) - 2074*pow(z,3) + 1538*pow(z,4) - 860*pow(z,5) + 107*pow(z,6))*log(z))/(pow(-1 + z,2)*pow(z,2)) + (36*(-35 + 280*z + 214*pow(z,2) - 280*pow(z,3) + 35*pow(z,4))*log(1 - z)*log(z))/pow(z,2) - 7560*(2*Li_2 + pow(log(1/(1 - z)),2) - pow(log(1 - z),2))*log(z) - (630*(-1 + z)*(1 + z)*(1 - 8*z + pow(z,2))*pow(log(z),2))/pow(z,2) - 2520*pow(log(z),3)))/(2016.*pow(C1,5))
//     ;

    
//  //------------------

//     H1_ext[1] = 
//     (-5*(1 + 12*pow(z,2) - 16*pow(z,3) + 3*pow(z,4) + 12*pow(z,2)*log(z)))/(32.*pow(C1,5)*pow(z,2)) - (pow(epsilon,2)*(-3727 + 31563*z - 62725*pow(z,2) + 2568*pow(Pi,2)*pow(z,2) + 76605*pow(z,3) - 7704*pow(Pi,2)*pow(z,3) - 79161*pow(z,4) + 7704*pow(Pi,2)*pow(z,4) + 54445*pow(z,5) - 2568*pow(Pi,2)*pow(z,5) - 19571*pow(z,6) + 2571*pow(z,7) + 30240*pow(z,2)*zeta_3 - 90720*pow(z,3)*zeta_3 + 90720*pow(z,4)*zeta_3 - 30240*pow(z,5)*zeta_3 + 30240*pow(-1 + z,3)*pow(z,2)*Li_3 + 1260*z*log(1/(1 - z)) - 12600*pow(z,2)*log(1/(1 - z)) + 6012*pow(z,3)*log(1/(1 - z)) + 30816*pow(z,4)*log(1/(1 - z)) - 36828*pow(z,5)*log(1/(1 - z)) + 12600*pow(z,6)*log(1/(1 - z)) - 1260*pow(z,7)*log(1/(1 - z)) - 630*pow(log(1/(1 - z)),2) + 1890*z*pow(log(1/(1 - z)),2) - 2034*pow(z,2)*pow(log(1/(1 - z)),2) + 11142*pow(z,3)*pow(log(1/(1 - z)),2) - 32562*pow(z,4)*pow(log(1/(1 - z)),2) + 36054*pow(z,5)*pow(log(1/(1 - z)),2) - 15750*pow(z,6)*pow(log(1/(1 - z)),2) + 1890*pow(z,7)*pow(log(1/(1 - z)),2) + 642*log(1 - z) - 666*z*log(1 - z) - 2970*pow(z,2)*log(1 - z) - 28014*pow(z,3)*log(1 - z) + 86670*pow(z,4)*log(1 - z) - 81126*pow(z,5)*log(1 - z) + 28650*pow(z,6)*log(1 - z) - 3186*pow(z,7)*log(1 - z) + 630*pow(log(1 - z),2) - 1890*z*pow(log(1 - z),2) + 2034*pow(z,2)*pow(log(1 - z),2) - 11142*pow(z,3)*pow(log(1 - z),2) + 32562*pow(z,4)*pow(log(1 - z),2) - 36054*pow(z,5)*pow(log(1 - z),2) + 15750*pow(z,6)*pow(log(1 - z),2) - 1890*pow(z,7)*pow(log(1 - z),2) - 840*log(z) + 11340*z*log(z) - 4194*pow(z,2)*log(z) + 8658*pow(z,3)*log(z) - 21672*pow(z,4)*log(z) + 24972*pow(z,5)*log(z) - 12270*pow(z,6)*log(z) + 1926*pow(z,7)*log(z) + 15120*pow(z,3)*log(1/(1 - z))*log(z) - 30240*pow(z,4)*log(1/(1 - z))*log(z) + 15120*pow(z,5)*log(1/(1 - z))*log(z) + 7560*pow(z,2)*pow(log(1/(1 - z)),2)*log(z) - 22680*pow(z,3)*pow(log(1/(1 - z)),2)*log(z) + 22680*pow(z,4)*pow(log(1/(1 - z)),2)*log(z) - 7560*pow(z,5)*pow(log(1/(1 - z)),2)*log(z) - 1260*log(1 - z)*log(z) + 3780*z*log(1 - z)*log(z) - 26604*pow(z,2)*log(1 - z)*log(z) + 105012*pow(z,3)*log(1 - z)*log(z) - 162972*pow(z,4)*log(1 - z)*log(z) + 109764*pow(z,5)*log(1 - z)*log(z) - 31500*pow(z,6)*log(1 - z)*log(z) + 3780*pow(z,7)*log(1 - z)*log(z) - 7560*pow(z,2)*pow(log(1 - z),2)*log(z) + 22680*pow(z,3)*pow(log(1 - z),2)*log(z) - 22680*pow(z,4)*pow(log(1 - z),2)*log(z) + 7560*pow(z,5)*pow(log(1 - z),2)*log(z) + 630*pow(log(z),2) - 1890*z*pow(log(z),2) + 9450*pow(z,2)*pow(log(z),2) - 33390*pow(z,3)*pow(log(z),2) + 54810*pow(z,4)*pow(log(z),2) - 43470*pow(z,5)*pow(log(z),2) + 15750*pow(z,6)*pow(log(z),2) - 1890*pow(z,7)*pow(log(z),2) + 2520*pow(z,2)*pow(log(z),3) - 7560*pow(z,3)*pow(log(z),3) + 7560*pow(z,4)*pow(log(z),3) - 2520*pow(z,5)*pow(log(z),3) + 36*pow(-1 + z,3)*Li_2*(35 + 848*pow(z,2) - 560*pow(z,3) + 105*pow(z,4) - 420*pow(z,2)*log(z))))/(2016.*pow(C1,5)*pow(-1 + z,3)*pow(z,2))
//     ;
    
//     H1_ext[1] *= dzdxi;

/*----------------------------------------*/
//Wrong Choice
/*----------------------------------------*/
//     H0_ext[0] = 
//     z + (pow(epsilon,2)*(177 - 4500*z - 6490*pow(z,2) + 420*pow(Pi,2)*pow(z,2) + 22522*pow(z,3) - 840*pow(Pi,2)*pow(z,3) - 12579*pow(z,4) + 420*pow(Pi,2)*pow(z,4) + 210*pow(z,5) - 2520*pow(-1 + z,2)*pow(z,2)*Li_2 + 6*pow(-1 + z,2)*(-35 + 280*z + 214*pow(z,2) - 280*pow(z,3) + 35*pow(z,4))*log(1 - z) + 210*log(z) - 2100*z*log(z) + 2286*pow(z,2)*log(z) + 2568*pow(z,3)*log(z) - 4854*pow(z,4)*log(z) + 2100*pow(z,5)*log(z) - 210*pow(z,6)*log(z) - 1260*pow(z,2)*pow(log(z),2) + 2520*pow(z,3)*pow(log(z),2) - 1260*pow(z,4)*pow(log(z),2)))/(630.*pow(-1 + z,2)*z)
//     ;
    
// //------------------
//     H1_ext[0] = 
//     1 + (pow(epsilon,2)*(-33 + 1569*z + 13204*pow(z,2) - 420*pow(Pi,2)*pow(z,2) - 41122*pow(z,3) + 1260*pow(Pi,2)*pow(z,3) + 42591*pow(z,4) - 1260*pow(Pi,2)*pow(z,4) - 15519*pow(z,5) + 420*pow(Pi,2)*pow(z,5) + 630*pow(z,6) - 2520*pow(-1 + z,3)*pow(z,2)*Li_2 + 6*pow(-1 + z,3)*(35 + 634*pow(z,2) - 560*pow(z,3) + 105*pow(z,4))*log(1 - z) + 210*log(z) - 630*z*log(z) + 4434*pow(z,2)*log(z) - 14982*pow(z,3)*log(z) + 22122*pow(z,4)*log(z) - 15774*pow(z,5)*log(z) + 5250*pow(z,6)*log(z) - 630*pow(z,7)*log(z) + 1260*pow(z,2)*pow(log(z),2) - 3780*pow(z,3)*pow(log(z),2) + 3780*pow(z,4)*pow(log(z),2) - 1260*pow(z,5)*pow(log(z),2)))/(630.*pow(-1 + z,3)*pow(z,2))
//     ;


//     H1_ext[0] *= dzdxi;

// //------------------

//     H0_ext[1] = (-5*(-1 + 8*z - 8*pow(z,3) + pow(z,4) + 12*pow(z,2)*log(z)))/(32.*pow(C1,5)*z) + (pow(epsilon,2)*(3871 - 105*pow(Pi,2) - 19501*z + 1050*pow(Pi,2)*z + 26454*pow(z,2) - 501*pow(Pi,2)*pow(z,2) + 12950*pow(z,3) - 2568*pow(Pi,2)*pow(z,3) - 47279*pow(z,4) + 3069*pow(Pi,2)*pow(z,4) + 26151*pow(z,5) - 1050*pow(Pi,2)*pow(z,5) - 2646*pow(z,6) + 105*pow(Pi,2)*pow(z,6) + 15120*pow(z,2)*zeta_3 - 30240*pow(z,3)*zeta_3 + 15120*pow(z,4)*zeta_3 - 15120*pow(-1 + z,2)*pow(z,2)*Li_3 - 321*log(1 - z) + 3210*z*log(1 - z) - 5457*pow(z,2)*log(1 - z) + 5457*pow(z,4)*log(1 - z) - 3210*pow(z,5)*log(1 - z) + 321*pow(z,6)*log(1 - z) - 210*log(z) + 10290*z*log(z) - 38577*pow(z,2)*log(z) + 1260*pow(Pi,2)*pow(z,2)*log(z) + 59442*pow(z,3)*log(z) - 2520*pow(Pi,2)*pow(z,3)*log(z) - 31224*pow(z,4)*log(z) + 1260*pow(Pi,2)*pow(z,4)*log(z) + 2580*pow(z,5)*log(z) - 321*pow(z,6)*log(z) + 630*log(1 - z)*log(z) - 6300*z*log(1 - z)*log(z) + 6858*pow(z,2)*log(1 - z)*log(z) + 7704*pow(z,3)*log(1 - z)*log(z) - 14562*pow(z,4)*log(1 - z)*log(z) + 6300*pow(z,5)*log(1 - z)*log(z) - 630*pow(z,6)*log(1 - z)*log(z) - 315*pow(log(z),2) + 3150*z*pow(log(z),2) - 5355*pow(z,2)*pow(log(z),2) + 5355*pow(z,4)*pow(log(z),2) - 3150*pow(z,5)*pow(log(z),2) + 315*pow(z,6)*pow(log(z),2) + 1260*pow(z,2)*pow(log(z),3) - 2520*pow(z,3)*pow(log(z),3) + 1260*pow(z,4)*pow(log(z),3) - 18*pow(-1 + z,2)*Li_2*(-35 + 280*z + 428*pow(z,2) - 280*pow(z,3) + 35*pow(z,4) - 420*pow(z,2)*log(z))))/(1008.*pow(C1,5)*pow(-1 + z,2)*z)
//     ;
    

    
//  //------------------

//     H1_ext[1] = 
//     (-5*(1 + 12*pow(z,2) - 16*pow(z,3) + 3*pow(z,4) + 12*pow(z,2)*log(z)))/(32.*pow(C1,5)*pow(z,2)) - (pow(epsilon,2)*(-4081 + 105*pow(Pi,2) + 22434*z - 315*pow(Pi,2)*z - 64625*pow(z,2) + 2859*pow(Pi,2)*pow(z,2) + 155830*pow(z,3) - 9417*pow(Pi,2)*pow(z,3) - 232503*pow(z,4) + 12987*pow(Pi,2)*pow(z,4) + 180230*pow(z,5) - 8529*pow(Pi,2)*pow(z,5) - 65223*pow(z,6) + 2625*pow(Pi,2)*pow(z,6) + 7938*pow(z,7) - 315*pow(Pi,2)*pow(z,7) + 15120*pow(z,2)*zeta_3 - 45360*pow(z,3)*zeta_3 + 45360*pow(z,4)*zeta_3 - 15120*pow(z,5)*zeta_3 + 15120*pow(-1 + z,3)*pow(z,2)*Li_3 + 321*log(1 - z) - 963*z*log(1 - z) + 4815*pow(z,2)*log(1 - z) - 17013*pow(z,3)*log(1 - z) + 27927*pow(z,4)*log(1 - z) - 22149*pow(z,5)*log(1 - z) + 8025*pow(z,6)*log(1 - z) - 963*pow(z,7)*log(1 - z) - 420*log(z) + 5670*z*log(z) - 28707*pow(z,2)*log(z) + 1260*pow(Pi,2)*pow(z,2)*log(z) + 84159*pow(z,3)*log(z) - 3780*pow(Pi,2)*pow(z,3)*log(z) - 90666*pow(z,4)*log(z) + 3780*pow(Pi,2)*pow(z,4)*log(z) + 39096*pow(z,5)*log(z) - 1260*pow(Pi,2)*pow(z,5)*log(z) - 6135*pow(z,6)*log(z) + 963*pow(z,7)*log(z) - 630*log(1 - z)*log(z) + 1890*z*log(1 - z)*log(z) - 13302*pow(z,2)*log(1 - z)*log(z) + 44946*pow(z,3)*log(1 - z)*log(z) - 66366*pow(z,4)*log(1 - z)*log(z) + 47322*pow(z,5)*log(1 - z)*log(z) - 15750*pow(z,6)*log(1 - z)*log(z) + 1890*pow(z,7)*log(1 - z)*log(z) + 315*pow(log(z),2) - 945*z*pow(log(z),2) + 4725*pow(z,2)*pow(log(z),2) - 16695*pow(z,3)*pow(log(z),2) + 27405*pow(z,4)*pow(log(z),2) - 21735*pow(z,5)*pow(log(z),2) + 7875*pow(z,6)*pow(log(z),2) - 945*pow(z,7)*pow(log(z),2) + 1260*pow(z,2)*pow(log(z),3) - 3780*pow(z,3)*pow(log(z),3) + 3780*pow(z,4)*pow(log(z),3) - 1260*pow(z,5)*pow(log(z),3) + 18*pow(-1 + z,3)*Li_2*(35 + 848*pow(z,2) - 560*pow(z,3) + 105*pow(z,4) - 420*pow(z,2)*log(z))))/(1008.*pow(C1,5)*pow(-1 + z,3)*pow(z,2))
//     ;
//     H1_ext[1] *= dzdxi;
    
}
//==========================================================
SERIES_ORIGIN_BULK::SERIES_ORIGIN_BULK(double l, double Omega, double nu0, double n, double brel, double gamma_0, 
    double fzeta_c, double H0_0, double W0_0):
l{l},
Omega{Omega},
nu0{nu0},
n{n},
brel{brel},
gamma_0{gamma_0},
fzeta_c{fzeta_c},
H0_0{H0_0},
W0_0{W0_0}
{

}

//-------------------------------------------------
SERIES_ORIGIN_BULK::~SERIES_ORIGIN_BULK(void)
{
}
//-------------------------------------------------
void SERIES_ORIGIN_BULK::operator()( double xi, const std::vector<double> yin, std::vector<double> &y )
{   
    
double H_0= yin[0];
double W_0= yin[1];
//-----------------------------
y[0] = H_0 + ((1 + n)*pow(xi,2)*((-3*pow(E,nu0)*fzeta_c*H0_0*Omega)/(3 + 2*l) + (fzeta_c*(1 + n)*W0_0*Omega*((1 + 3*brel)*pow(E,nu0)*l - 3*pow(Omega,2)))/(l*(3 + 2*l)) - (2*pow(E,nu0/2.)*H_0*gamma_0*(pow(E,nu0)*(1 + l)*(3 + 3*pow(brel,2)*(9 + l)*gamma_0 + brel*(3 + (15 - 3*l - 2*pow(l,2))*gamma_0)) + 3*brel*(9 + l)*gamma_0*pow(Omega,2)))/((1 + l)*(6 + 4*l)) - ((1 + brel)*pow(E,nu0/2.)*W_0*gamma_0*((1 + 3*brel)*pow(E,nu0)*l*(1 + l)*(-1 + n*(-1 + gamma_0)) - 3*(1 + n)*(-1 + 9*brel*gamma_0 + l*(-1 + brel*gamma_0))*pow(Omega,2)))/(l*(1 + l)*(3 + 2*l))))/(6.*pow(E,(3*nu0)/2.)*pow(gamma_0,2));
//----------------------
y[1] = W_0 + (pow(xi,2)*((-3*fzeta_c*H0_0*Omega*((1 + 3*brel)*pow(E,nu0)*l*(1 + l) + 3*(2 + l)*pow(Omega,2)))/((1 + brel)*pow(E,nu0/2.)) + 3*H_0*gamma_0*(((1 + 3*brel)*pow(E,nu0)*l*(1 + l)*(-1 + n*(-1 + gamma_0)))/(1 + n) + 3*(-2 - 6*brel*gamma_0 + l*(-1 + brel*gamma_0))*pow(Omega,2)) + (fzeta_c*(1 + n)*W0_0*Omega*(pow(E,2*nu0)*(1 + l)*pow(l + 3*brel*l,2) + 3*(1 + 3*brel)*pow(E,nu0)*l*pow(Omega,2) - 9*(2 + l)*pow(Omega,4)))/((1 + brel)*pow(E,(3*nu0)/2.)*l) + (W_0*gamma_0*(-(pow(E,2*nu0)*(1 + l)*pow(l + 3*brel*l,2)*(-1 + n*(-1 + gamma_0))) + 3*pow(E,nu0)*l*(1 + n + n*gamma_0 + l*n*gamma_0 + 6*pow(brel,2)*(-1 + l)*(1 + n)*gamma_0 + brel*(3 + 2*(-4 + 3*l + pow(l,2))*gamma_0 + n*(3 + (-5 + 9*l + 2*pow(l,2))*gamma_0)))*pow(Omega,2) - 9*(2 + l)*(1 + n)*pow(Omega,4)))/(pow(E,nu0)*l)))/(18.*(3 + 2*l)*pow(gamma_0,2)*pow(Omega,2));
//----------------------
y[2] = -(W_0/l) + (pow(xi,2)*((3*fzeta_c*H0_0*Omega*((1 + 3*brel)*pow(E,nu0)*(3 + l) + 3*pow(Omega,2)))/((1 + brel)*pow(E,nu0/2.)) - (3*H_0*gamma_0*((1 + 3*brel)*pow(E,nu0)*(3 + 4*l + pow(l,2))*(-1 + n*(-1 + gamma_0)) + 3*(1 + n)*(-1 + 9*brel*gamma_0 + l*(-1 + brel*gamma_0))*pow(Omega,2)))/((1 + l)*(1 + n)) - (fzeta_c*(1 + n)*W0_0*Omega*(pow(1 + 3*brel,2)*pow(E,2*nu0)*l*(3 + l) - 9*(1 + 3*brel)*pow(E,nu0)*pow(Omega,2) - 9*pow(Omega,4)))/((1 + brel)*pow(E,(3*nu0)/2.)*l) + (W_0*gamma_0*(pow(1 + 3*brel,2)*pow(E,2*nu0)*l*(3 + 4*l + pow(l,2))*(-1 + n*(-1 + gamma_0)) - 3*pow(E,nu0)*(6*pow(brel,2)*(-3 + 2*l + pow(l,2))*(1 + n)*gamma_0 + (1 + l)*(-3 + n*(-3 + (3 + l)*gamma_0)) + brel*(2*pow(l,3)*(1 + n)*gamma_0 + pow(l,2)*(8 + 11*n)*gamma_0 + 3*l*(-3 - 3*n + 4*n*gamma_0) - 3*(3 + 3*n + 10*gamma_0 + 7*n*gamma_0)))*pow(Omega,2) + 9*(1 + l)*(1 + n)*pow(Omega,4)))/(pow(E,nu0)*l*(1 + l))))/(18.*(3 + 2*l)*pow(gamma_0,2)*pow(Omega,2));
//----------------------
y[3] = ((1 + n)*xi*((-3*pow(E,nu0)*fzeta_c*H0_0*Omega)/(3 + 2*l) + (fzeta_c*(1 + n)*W0_0*Omega*((1 + 3*brel)*pow(E,nu0)*l - 3*pow(Omega,2)))/(l*(3 + 2*l)) - (2*pow(E,nu0/2.)*H_0*gamma_0*(pow(E,nu0)*(1 + l)*(3 + 3*pow(brel,2)*(9 + l)*gamma_0 + brel*(3 + (15 - 3*l - 2*pow(l,2))*gamma_0)) + 3*brel*(9 + l)*gamma_0*pow(Omega,2)))/((1 + l)*(6 + 4*l)) - ((1 + brel)*pow(E,nu0/2.)*W_0*gamma_0*((1 + 3*brel)*pow(E,nu0)*l*(1 + l)*(-1 + n*(-1 + gamma_0)) - 3*(1 + n)*(-1 + 9*brel*gamma_0 + l*(-1 + brel*gamma_0))*pow(Omega,2)))/(l*(1 + l)*(3 + 2*l))))/(3.*pow(E,(3*nu0)/2.)*pow(gamma_0,2));
//----------------------
    
}
//-------------------------------------------------