#ifndef _LE_SOLVER_HPP_
#define _LE_SOLVER_HPP_

#include <vector>
// #include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>


class LE_SOLVER{
    typedef std::vector<double> state_type;
    // boost::numeric::odeint::runge_kutta4< state_type > bs;
    boost::numeric::odeint::bulirsch_stoer< state_type > bs{ 1e-12 , 1e-12 };
public:
    LE_SOLVER(double n, double b);
    ~LE_SOLVER(void);
    double n;
    double b;
    double nu0;
//------------------------------------------------
    void solve_adaptive(const double ximin, const int N, double &xiguess, double &theta_f, double &mu_f, double &nu_f);

    void solve_and_interp(
    int N, double xil, double xir, double &dxi,
    double theta_f, double mu_f, double nu_f,
    std::vector<double> &xivals,
    std::vector<double> &theta,
    std::vector<double> &lambda,
    std::vector<double> &nu);

    void step_controlled(state_type &y, double &xi0, double &dxi, double &dximax);
//-------------------------------------------------
private:
//---------------------------------------------------------------
    double dthetadxi(double xi, double theta, double mu, double nu);
    double dmudxi(double xi, double theta, double mu, double nu);
    double dnudxi(double xi, double theta, double mu, double nu);

    void series_origin(double xi, double &theta, double &mu, double &nu);
    
    void LE_rhs( const state_type &y , state_type &dydt , double xi ); 
    void step(state_type &y, double &xi, double &dxi);
    
};




#endif