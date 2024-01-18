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
//-------------------------------------------------
#include "lane-emden-solver.hpp"
#include "outputfiles.hpp"
//-------------------------------------------------

LE_SOLVER::LE_SOLVER(double n, double b):
n{n},
b{b},
nu0{0.}
{
}
//-------------------------------------------------
LE_SOLVER::~LE_SOLVER(void)
{

}
//-------------------------------------------------
double LE_SOLVER::dthetadxi(double xi, double theta, double mu, double nu)
{   
    return 
    -(((1 + b*theta)*(b*pow(xi,3)*pow(abs(theta),1 + n) + mu))/(xi*(xi - 2*b*(1 + n)*mu)))
    ;
}
double LE_SOLVER::dmudxi(double xi, double theta, double mu, double nu)
{
    return (xi*xi)*pow(abs(theta),n);
}
double LE_SOLVER::dnudxi(double xi, double theta, double mu, double nu)
{
    return 
    (-2*b*(1 + n)*(b*pow(xi,3)*pow(abs(theta),1 + n) + mu))/(xi*(-xi + 2*b*(1 + n)*mu))
    ;
}
void LE_SOLVER::series_origin(double xi, double &theta, double &mu, double &nu)
{   
    double nu0 = 1.;

    theta = 1 - ((1 + b)*(1 + 3*b)*pow(xi,2))/6. + ((1 + b)*(1 + 3*b)*(30*pow(b,2) + 3*n - 2*b*n + 15*pow(b,2)*n)*pow(xi,4))/360. - ((1 + b)*(1 + 3*b)*(210*pow(b,2) + 1890*pow(b,4) - 15*n - 180*b*n + 398*pow(b,2)*n - 312*pow(b,3)*n + 2205*pow(b,4)*n + 24*pow(n,2) - 90*b*pow(n,2) + 130*pow(b,2)*pow(n,2) - 246*pow(b,3)*pow(n,2) + 630*pow(b,4)*pow(n,2))*pow(xi,6))/45360.
    ;

    mu = pow(xi,3)/3. - ((1 + b)*(1 + 3*b)*n*pow(xi,5))/30. + ((1 + b)*(1 + 3*b)*n*(-5 - 20*b + 15*pow(b,2) + 8*n + 18*b*n + 30*pow(b,2)*n)*pow(xi,7))/2520.
    ;

    nu = 
    nu0 + (b*(1 + 3*b)*(1 + n)*pow(xi,2))/3. - (b*(1 + 3*b)*(1 + n)*(-5*b + 15*pow(b,2) + 3*n - 2*b*n + 15*pow(b,2)*n)*pow(xi,4))/180. + (b*(1 + 3*b)*(1 + n)*(280*pow(b,2) - 210*pow(b,3) + 630*pow(b,4) - 15*n - 243*b*n + 251*pow(b,2)*n - 501*pow(b,3)*n + 1260*pow(b,4)*n + 24*pow(n,2) - 90*b*pow(n,2) + 130*pow(b,2)*pow(n,2) - 246*pow(b,3)*pow(n,2) + 630*pow(b,4)*pow(n,2))*pow(xi,6))/22680.
    ;
}
//-------------------------------------------------
void LE_SOLVER::LE_rhs( const state_type &y , state_type &dydt , double xi )
{
    double theta = y[0], mu = y[1], nu = y[2];

    dydt[0] = dthetadxi(xi, theta, mu, nu);
    dydt[1] = dmudxi(xi,theta,mu, nu);
    dydt[2] = dnudxi(xi, theta, mu, nu);
}
// //-------------------------------------------------
void LE_SOLVER::step(state_type &y, double &xi, double &dxi)
{   
    // cout<<"xi = "<<xi<<endl;
    // cout<<"y[0] = "<<y[0]<<"; y[1] = "<<y[1]<<endl;
    double tol = 1e-10;
    state_type yin = y;
    double xiin = xi;
    bs.try_step([this] (const state_type &y , state_type &dydt , double xi ) 
    { return LE_rhs (y, dydt,xi); }, y, xi,dxi);
    
    if(y[0]>tol)
    {

    }
    else{
        if(y[0]>0. || abs(y[0])<tol)
        {
        }
        else{
            double factor = 2.;
            int max_steps = 100000;
            double y0 = y[0];
            int counter = 0;
            double dx = dxi;
            while(counter<max_steps && y0<0. && abs(y0)>tol )
            {   
                // cout<<"beginning step"<<endl;
                // cout<<"y0 = "<<y0<<endl;
                // cout<<"dx = "<<dx<<endl;
                // cout<<"xi = "<<xi<<endl;
                xi = xiin;
                y = yin;
                dx = dxi/(factor*(counter+1.));
                bs.try_step([this] (const state_type &y , state_type &dydt , double xi ) 
                 { return LE_rhs (y, dydt,xi); }, y, xi,dx);
                y0 = y[0];
                // cout<<"y0 = "<<y0<<endl;
                // cout<<"dx = "<<dx<<endl;
                // cout<<"xi = "<<xi<<endl;
                counter +=1;
            }
            // cout<<"end step"<<endl;
            dxi = dx;

        }
    }
}
// //-------------------------------------------------------------
void LE_SOLVER::solve_adaptive(const double ximin, const int N, double &xiguess, double &theta_f, double &mu_f, double &nu_f)
{
    double dxi = (4.-ximin)/N;
    double theta0 =0., mu0 = 0., nu0 = 0.;
    series_origin(ximin, theta0,mu0, nu0);
    vector<double> yvals = {theta0,mu0, nu0};
    int Nmax = 0;
    double xi = ximin;
    double thetaf=0., muf = 0., nuf= 0.;

    for(int i = 1; i<N; ++i)
    {   
        thetaf = yvals[0]; muf = yvals[1]; nuf = yvals[2];
        step(yvals, xi, dxi ); 
        if(yvals[0]<0.)
        {   
            Nmax = i-1;
            // cout<<"theta = "<<yvals[0]<<"; theta_im1 = "<<theta[i-1]<<endl;
            // cout<<"xi1 = "<<xi<<"; mu1 = "<<mu[i-1]<<endl;
            break;
        }
        
    }
    xiguess = xi;
    theta_f = thetaf;
    mu_f = muf;
    nu_f = nuf;
    cout<<"xi1 = "<<xi<<endl;
    cout<<"thetai = "<<thetaf<<"; mu1 = "<<muf<<"; nuf = "<<nu_f<<endl; 
    cout<<"2M/R = "<<2*b*(n+1)*mu_f/xi<<endl;
    cout<<"Nmax = "<<Nmax <<endl;
    cout<<"-----------"<<endl;
    // cout<<"updating xi"<<endl;
    // double a = muf/(xi*xi);
    // double xiupdate = thetaf/(a) + xi;
    // xiguess = xiupdate;
    // cout<<"diff = "<<xiguess-xi<<endl;

    
}
// //-------------------------------------------------------------
// //-------------------------------------------------------------
void LE_SOLVER::step_controlled(state_type &yvals, double &xi0, double &dxi, double &dximax)
{   
    if(fabs(dxi)>fabs(dximax))
    {
        dxi = dximax/2.;
    } 
    if(dximax<0){ dximax = -dximax; }
    double xif = xi0 + dximax;
    assert(fabs(dxi) < fabs(dximax));
    assert(xif>xi0);
    assert(xi0>0.);

    if(dxi<0){ dxi = -dxi; }
    double xi = xi0;
    int counter = 0;
    while(counter == 0)
    {       
        vector<double> y = yvals;
        double xi_0 = xi;
        double dxi0 = dxi;
        step(y, xi_0, dxi0);
        
        if(xi_0<xif)
        {
            yvals = y;
            xi = xi_0;
            dxi = dxi0;
        }
        else{
            counter = 1;
        }
        
        
    }
        dxi = xif - xi;
        step(yvals, xi, dxi);
        xi0 = xi;
    }
// //---------------------------------------------------------
void LE_SOLVER::solve_and_interp(
    int N, double xil, double xir, double &dxi,
    double theta_f, double mu_f, double nu_f,
    std::vector<double> &xivals,
    std::vector<double> &theta,
    std::vector<double> &lambda,
    std::vector<double> &nu)
{   
    {
        int N1 = theta.size();
        assert(N1==N);
        N1 = lambda.size();
        assert(N1==N);
        N1 = nu.size();
        assert(N1==N);
    }
    std::vector<double> mu(N);
    
    double dxil = (xir - xil)/(N-1);

    // cout<<"dxil = "<<dxil<<endl;

    double dxil_start = dxil/2.;

    double xil_start = xil;

    {
        int i=0;

        series_origin(xil_start, theta[i], mu[i], nu[i]);

        lambda[i] = -log(1.0 - 2.0*b*mu[i]*(1+n)/xil_start);
        

        xivals[i] = xil_start;

        // cout<<"xivals = "<<xivals[i]<<endl;
        // cout<<"theta = "<<theta[i]<<endl;
        // cout<<"mu = "<<mu[i]<<endl;

    }

    for(int i=1; i< N-1; ++i)
    {
        vector<double> yvals(3);
        yvals[0] = theta[i-1];
        yvals[1] = mu[i-1];
        yvals[2] = nu[i-1];

        // cout<<"xivals = "<<xivals[i]<<endl;
        // cout<<"theta = "<<theta[i]<<endl;
        // cout<<"mu = "<<mu[i]<<endl;


        step_controlled(yvals, xil_start, dxil_start, dxil);

        theta[i] = yvals[0];
        mu[i] = yvals[1];
        nu[i] = yvals[2];

        lambda[i] = -log(1.0 - 2.0*b*mu[i]*(1+n)/xil_start);

        xivals[i] = xil_start;

        // cout<<"i = "<<i<<endl;
        // cout<<"xivals = "<<xivals[i]<<endl;
        // cout<<"theta = "<<theta[i]<<endl;
        // cout<<"mu = "<<mu[i]<<endl;
        // cout<<"diff xi = "<<xivals[i]-xivals[i-1]<<endl;
        // cout<<"---------"<<endl;


    }
    {
        int i = N-1;

        theta[i] = theta_f;
        mu[i] = mu_f;
        nu[i] = nu_f;

        lambda[i] = -log(1.0 - 2.0*b*mu[i]*(1+n)/xir);

        xivals[i] = xir;

        // cout<<"i = "<<i<<endl;
        // cout<<"xivals = "<<xivals[i]<<endl;
        // cout<<"theta = "<<theta[i]<<endl;
        // cout<<"mu = "<<mu[i]<<endl;
        // cout<<"---------"<<endl;

    }
    dxi = dxil;
    cout<<"dxi = "<<dxi<<endl;

    double shift = nu_f - log(1 - (2.*b*(1+n)*mu_f)/(xir));

    shift_vec(shift, nu);

    nu0 = nu[0];
}