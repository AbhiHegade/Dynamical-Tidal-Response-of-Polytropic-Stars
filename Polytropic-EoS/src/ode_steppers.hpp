// #ifndef _ODE_STEPPERS_HPP_
// #define _ODE_STEPPERS_HPP_

// #include <vector>
// #include <string>
// #include <functional>


// void step(std::function <void(const std::vector<double>& , std::vector<double>& , double  )> rhs, 
// std::vector<double> &y, double &r, double &dr);

// void step_controlled(std::function <void(const std::vector<double>& , std::vector<double>& , double )> rhs, 
// std::string side, std::vector<double> &yvals, double &r0, double &dr, double &drmax);

// void solve_adaptive(
//     std::function <void(const std::vector<double>& , std::vector<double>& , double )> rhs ,
//     std::function <void(double , const std::vector<double> , std::vector<double>& )> series_origin,  
//     std::function <void(double , const std::vector<double> , std::vector<double>& )> series_surf,      
//     std::string side, int get_series,  double r0, double rf, double R_star,
//     const std::vector<double> ics,
//     std::vector<double> &series_sol, 
//     std::vector<double> &yvals,
//     double &dr);



// #endif

#ifndef _ODE_STEPPERS_HPP_
#define _ODE_STEPPERS_HPP_

#include <vector>
#include <string>
#include <functional>


// void step(std::function <void(const std::vector<double>& , std::vector<double>& , double  )> rhs, 
// std::vector<double> &y, double &r, double &dr);

void step_controlled(std::function <void(const std::vector<double>& , std::vector<double>& , double )> rhs, 
std::string side, std::vector<double> &yvals, double &r0, double &dr, double &drmax);

void solve_adaptive(
    std::function <void(const std::vector<double>& , std::vector<double>& , double )> rhs ,
    std::function <void(double , const std::vector<double> , std::vector<double>& )> series_origin,  
    std::function <void(double , const std::vector<double> , std::vector<double>& )> series_surf,      
    std::string side, int get_series,  double r0, double rf, double R_star,
    const std::vector<double> ics,
    std::vector<double> &series_sol, 
    std::vector<double> &yvals,
    double &dr);


struct push_back_state_and_time
{
    std::vector< std::vector<double> > &m_states;

    push_back_state_and_time( std::vector< std::vector<double> > &states);

    // : m_states( states ) , m_times( times ) { }

    void operator()( const std::vector<double> &x , double t );
    
};

void solve_adaptive_and_write(std::function <void(const std::vector<double>& , std::vector<double>& , double )> rhs,
    std::vector<double> &rl,
    std::vector<double> icsl,
    std::vector<std::vector<double>> &sol,
    std::string side
);





#endif