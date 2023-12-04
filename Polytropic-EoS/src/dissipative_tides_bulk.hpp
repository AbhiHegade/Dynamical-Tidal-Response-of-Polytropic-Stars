#ifndef _DISS_BULK_TIDES_HPP_
#define _DISS_BULK_TIDES_HPP_

#include <vector>
#include <string>
#include <functional>

#include "alpha-beta-funcs.hpp"
#include "local_expansions.hpp"

class DISS_BULK_TIDE{
public:
    DISS_BULK_TIDE(double l, double Omega, 
    double n, double brel,
    double nu0,
    double gamma_0,
    double mu1,
    double xi1,
    std::function <double(double)> thetafunc, 
    std::function <double(double)> lambdafunc, 
    std::function <double(double)> nufunc,
    std::function <double(double, double)> fzetafunc,
    std::function <double(double)> H00, 
    std::function <double(double)> W0, 
    std::function <double(double)> V0,
    double H0_0, double W0_0
        );
    ~DISS_BULK_TIDE(void);
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

    std::function <double(double, double)> fzetafunc;
    std::function <double(double)> H00func;
    std::function <double(double)> W0func;
    std::function <double(double)> V0func;

    double H0_0; double W0_0;

    RHS_PF_BULK rhspart;
    SERIES_ORIGIN_BULK series_origin_part;
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
// //=========================================================================
// void series_surf_0(
//         double r,
//         const std::vector<double> yin, std::vector<double> &y );
//==========================================================================
// void step_0(
//      state_type &y, double &r, double &dr);
//==========================================================================

};



































#endif
