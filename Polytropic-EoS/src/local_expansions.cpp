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
    z + (pow(epsilon,2)*(1197 - 22890*z + 37465*pow(z,2) + 1680*pow(Pi,2)*pow(z,2) - 20136*pow(z,3) - 3360*pow(Pi,2)*pow(z,3) - 3517*pow(z,4) + 1680*pow(Pi,2)*pow(z,4) + 5730*pow(z,5) - 489*pow(z,6) - 10080*pow(-1 + z,2)*pow(z,2)*Li_2 + 840*log(2) - 8400*z*log(2) + 9144*pow(z,2)*log(2) + 10272*pow(z,3)*log(2) - 19416*pow(z,4)*log(2) + 8400*pow(z,5)*log(2) - 840*pow(z,6)*log(2) + 24*pow(-1 + z,2)*(-35 + 280*z + 214*pow(z,2) - 280*pow(z,3) + 35*pow(z,4))*log(1 - z) + 840*log(z) - 8400*z*log(z) + 3276*pow(z,2)*log(z) + 22008*pow(z,3)*log(z) - 25284*pow(z,4)*log(z) + 8400*pow(z,5)*log(z) - 840*pow(z,6)*log(z) - 10080*pow(z,2)*log(2)*log(z) + 20160*pow(z,3)*log(2)*log(z) - 10080*pow(z,4)*log(2)*log(z) - 5040*pow(z,2)*pow(log(z),2) + 10080*pow(z,3)*pow(log(z),2) - 5040*pow(z,4)*pow(log(z),2)))/(2520.*pow(-1 + z,2)*z)
    ;

//------------------
    H1_ext[0] = 
    1 + (pow(epsilon,2)*(357 + 4809*z + 5039*pow(z,2) - 1680*pow(Pi,2)*pow(z,2) - 25069*pow(z,3) + 5040*pow(Pi,2)*pow(z,3) + 47571*pow(z,4) - 5040*pow(Pi,2)*pow(z,4) - 40705*pow(z,5) + 1680*pow(Pi,2)*pow(z,5) + 14745*pow(z,6) - 1467*pow(z,7) - 10080*pow(-1 + z,3)*pow(z,2)*Li_2 + 840*log(2) - 2520*z*log(2) + 17736*pow(z,2)*log(2) - 59928*pow(z,3)*log(2) + 88488*pow(z,4)*log(2) - 63096*pow(z,5)*log(2) + 21000*pow(z,6)*log(2) - 2520*pow(z,7)*log(2) + 24*pow(-1 + z,3)*(35 + 634*pow(z,2) - 560*pow(z,3) + 105*pow(z,4))*log(1 - z) + 840*log(z) - 2520*z*log(z) + 23604*pow(z,2)*log(z) - 77532*pow(z,3)*log(z) + 106092*pow(z,4)*log(z) - 68964*pow(z,5)*log(z) + 21000*pow(z,6)*log(z) - 2520*pow(z,7)*log(z) + 10080*pow(z,2)*log(2)*log(z) - 30240*pow(z,3)*log(2)*log(z) + 30240*pow(z,4)*log(2)*log(z) - 10080*pow(z,5)*log(2)*log(z) + 5040*pow(z,2)*pow(log(z),2) - 15120*pow(z,3)*pow(log(z),2) + 15120*pow(z,4)*pow(log(z),2) - 5040*pow(z,5)*pow(log(z),2)))/(2520.*pow(-1 + z,3)*pow(z,2))
    ;


    H1_ext[0] *= dzdxi;

//------------------

    H0_ext[1] = (-5*(-1 + 8*z - 8*pow(z,3) + pow(z,4) + 12*pow(z,2)*log(z)))/(32.*pow(C1,5)*z) + (pow(epsilon,2)*(125509 - 14700*pow(Pi,2) + 1434170*z + 147000*pow(Pi,2)*z - 3375767*pow(z,2) - 70140*pow(Pi,2)*pow(z,2) + 1813000*pow(z,3) - 359520*pow(Pi,2)*pow(z,3) + 460267*pow(z,4) + 429660*pow(Pi,2)*pow(z,4) - 503170*pow(z,5) - 147000*pow(Pi,2)*pow(z,5) + 45991*pow(z,6) + 14700*pow(Pi,2)*pow(z,6) + 2116800*pow(z,2)*zeta_3 - 4233600*pow(z,3)*zeta_3 + 2116800*pow(z,4)*zeta_3 - 2116800*pow(-1 + z,2)*pow(z,2)*Li_3 + 44940*log(2) - 449400*z*log(2) + 763980*pow(z,2)*log(2) - 763980*pow(z,4)*log(2) + 449400*pow(z,5)*log(2) - 44940*pow(z,6)*log(2) - 44940*log(1 - z) + 449400*z*log(1 - z) - 763980*pow(z,2)*log(1 - z) + 763980*pow(z,4)*log(1 - z) - 449400*pow(z,5)*log(1 - z) + 44940*pow(z,6)*log(1 - z) - 29400*log(z) + 1440600*z*log(z) - 403608*pow(z,2)*log(z) + 176400*pow(Pi,2)*pow(z,2)*log(z) - 1672464*pow(z,3)*log(z) - 352800*pow(Pi,2)*pow(z,3)*log(z) + 625812*pow(z,4)*log(z) + 176400*pow(Pi,2)*pow(z,4)*log(z) + 361200*pow(z,5)*log(z) - 44940*pow(z,6)*log(z) - 539280*pow(z,2)*log(2)*log(z) + 1078560*pow(z,3)*log(2)*log(z) - 539280*pow(z,4)*log(2)*log(z) + 88200*log(1 - z)*log(z) - 882000*z*log(1 - z)*log(z) + 960120*pow(z,2)*log(1 - z)*log(z) + 1078560*pow(z,3)*log(1 - z)*log(z) - 2038680*pow(z,4)*log(1 - z)*log(z) + 882000*pow(z,5)*log(1 - z)*log(z) - 88200*pow(z,6)*log(1 - z)*log(z) - 44100*pow(log(z),2) + 441000*z*pow(log(z),2) - 749700*pow(z,2)*pow(log(z),2) + 749700*pow(z,4)*pow(log(z),2) - 441000*pow(z,5)*pow(log(z),2) + 44100*pow(z,6)*pow(log(z),2) + 176400*pow(z,2)*pow(log(z),3) - 352800*pow(z,3)*pow(log(z),3) + 176400*pow(z,4)*pow(log(z),3) - 2520*pow(-1 + z,2)*Li_2*(-35 + 280*z + 428*pow(z,2) - 280*pow(z,3) + 35*pow(z,4) - 420*pow(z,2)*log(z))))/(141120.*pow(C1,5)*pow(-1 + z,2)*z)
    ;

    
 //------------------

    H1_ext[1] = 
    (-5*(1 + 12*pow(z,2) - 16*pow(z,3) + 3*pow(z,4) + 12*pow(z,2)*log(z)))/(32.*pow(C1,5)*pow(z,2)) - (pow(epsilon,2)*(-154909 + 14700*pow(Pi,2) + 1891467*z - 44100*pow(Pi,2)*z - 2801035*pow(z,2) + 400260*pow(Pi,2)*pow(z,2) - 254643*pow(z,3) - 1318380*pow(Pi,2)*pow(z,3) + 3679077*pow(z,4) + 1818180*pow(Pi,2)*pow(z,4) - 3501539*pow(z,5) - 1194060*pow(Pi,2)*pow(z,5) + 1279555*pow(z,6) + 367500*pow(Pi,2)*pow(z,6) - 137973*pow(z,7) - 44100*pow(Pi,2)*pow(z,7) + 2116800*pow(z,2)*zeta_3 - 6350400*pow(z,3)*zeta_3 + 6350400*pow(z,4)*zeta_3 - 2116800*pow(z,5)*zeta_3 + 2116800*pow(-1 + z,3)*pow(z,2)*Li_3 - 44940*log(2) + 134820*z*log(2) - 674100*pow(z,2)*log(2) + 2381820*pow(z,3)*log(2) - 3909780*pow(z,4)*log(2) + 3100860*pow(z,5)*log(2) - 1123500*pow(z,6)*log(2) + 134820*pow(z,7)*log(2) + 44940*log(1 - z) - 134820*z*log(1 - z) + 674100*pow(z,2)*log(1 - z) - 2381820*pow(z,3)*log(1 - z) + 3909780*pow(z,4)*log(1 - z) - 3100860*pow(z,5)*log(1 - z) + 1123500*pow(z,6)*log(1 - z) - 134820*pow(z,7)*log(1 - z) - 58800*log(z) + 793800*z*log(z) + 978192*pow(z,2)*log(z) + 176400*pow(Pi,2)*pow(z,2)*log(z) - 3209256*pow(z,3)*log(z) - 529200*pow(Pi,2)*pow(z,3)*log(z) + 2298276*pow(z,4)*log(z) + 529200*pow(Pi,2)*pow(z,4)*log(z) + 476268*pow(z,5)*log(z) - 176400*pow(Pi,2)*pow(z,5)*log(z) - 858900*pow(z,6)*log(z) + 134820*pow(z,7)*log(z) - 539280*pow(z,2)*log(2)*log(z) + 1617840*pow(z,3)*log(2)*log(z) - 1617840*pow(z,4)*log(2)*log(z) + 539280*pow(z,5)*log(2)*log(z) - 88200*log(1 - z)*log(z) + 264600*z*log(1 - z)*log(z) - 1862280*pow(z,2)*log(1 - z)*log(z) + 6292440*pow(z,3)*log(1 - z)*log(z) - 9291240*pow(z,4)*log(1 - z)*log(z) + 6625080*pow(z,5)*log(1 - z)*log(z) - 2205000*pow(z,6)*log(1 - z)*log(z) + 264600*pow(z,7)*log(1 - z)*log(z) + 44100*pow(log(z),2) - 132300*z*pow(log(z),2) + 661500*pow(z,2)*pow(log(z),2) - 2337300*pow(z,3)*pow(log(z),2) + 3836700*pow(z,4)*pow(log(z),2) - 3042900*pow(z,5)*pow(log(z),2) + 1102500*pow(z,6)*pow(log(z),2) - 132300*pow(z,7)*pow(log(z),2) + 176400*pow(z,2)*pow(log(z),3) - 529200*pow(z,3)*pow(log(z),3) + 529200*pow(z,4)*pow(log(z),3) - 176400*pow(z,5)*pow(log(z),3) + 2520*pow(-1 + z,3)*Li_2*(35 + 848*pow(z,2) - 560*pow(z,3) + 105*pow(z,4) - 420*pow(z,2)*log(z))))/(141120.*pow(C1,5)*pow(-1 + z,3)*pow(z,2))
    ;

    H1_ext[1] *= dzdxi;
    
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