#include <functional>
using std::function;
#include<string>
using std::string;
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
#include "ode_steppers.hpp"
#include "outputfiles.hpp"
//-------------------------------------------------
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
boost::numeric::odeint::bulirsch_stoer< std::vector<double>  > bs{ 1e-8 , 1e-8 };
//-------------------------------------------------
void step(std::function <void(const std::vector<double>& , std::vector<double>& , double  )> rhs, 
std::vector<double> &y, double &r, double &dr)
{
    bs.try_step(rhs, y, r,dr);
}

void step_controlled(std::function <void(const std::vector<double>& , std::vector<double>& , double )> rhs, 
std::string side, std::vector<double> &yvals, double &r0, double &dr, double &drmax)
{
    if(fabs(dr)>fabs(drmax))
    {
        dr = drmax/2.;
    }
    if(side == "left")
    {   
        if(drmax<0){ drmax = -drmax; }
        double rf = r0 + drmax;
        assert(fabs(dr) < fabs(drmax));
        assert(rf>r0);
        assert(r0>0.);
        if(dr<0){ dr = -dr; }
        double r = r0;
        int counter = 0;
        while(counter == 0)
        {       
            vector<double> y = yvals;
            double r_0 = r;
            double dr0 = dr;
           
            step(rhs, y, r_0, dr0);
                
            if(r_0<rf)
            {

                yvals = y;
                r = r_0;
                dr = dr0;
            }
            else{
                counter = 1;
            }
            
            
        }
            dr = rf - r;
            step(rhs, yvals, r, dr);
            r0 = r;
    }
    else if(side == "right"){
        if(drmax>0){ drmax = -drmax; }
        double rf = r0 + drmax;
        assert(fabs(dr) < fabs(drmax));
        assert(r0>rf);
        if(dr>0){ dr = -dr; }
        
        double r = r0;
        int counter = 0;
        while(counter ==0)
        {   
            vector<double> y = yvals;
            double r_0 = r;
            double dr0 = dr;
            step(rhs, y, r_0, dr0);
            if(r_0>rf)
            {
                yvals = y;
                r = r_0;
                dr = dr0;
            }
            else{
                counter = 1;
            }
            
        }
            dr = rf - r;
            step(rhs, yvals, r, dr);
            r0 = r;

    }
}

void solve_adaptive(
    std::function <void(const std::vector<double>& , std::vector<double>& , double )> rhs ,
    std::function <void(double , const std::vector<double> , std::vector<double>& )> series_origin,  
    std::function <void(double , const std::vector<double> , std::vector<double>& )> series_surf,      
    std::string side, int get_series,  double r0, double rf, double R_star,
    const std::vector<double> ics,
    std::vector<double> &series_sol, 
    std::vector<double> &yvals,
    double &dr)
{
    assert(yvals.size() == series_sol.size());
    
    if(side == "left")
    {   
        assert(rf>r0);
        assert(r0>0.);
        if(dr<0){ dr = -dr; }
        if(get_series)
        {   
            assert(ics.size()==2);
            series_origin(r0, ics, series_sol);

            // cout<<"get_series_true"<<endl;
            // cout<<side<<endl;
            

        }
        
        yvals = series_sol;
        // cout<<"part = "<<sol_type<<endl;
        // cout<<"series at origin"<<endl;
        // out_put_vec(yvals);
        // cout<<"---"<<endl;
        double r = r0;
        int counter = 0;
        while(counter == 0)
        {       
            vector<double> y = yvals;
            double r_0 = r;
            double dr0 = dr;

            // out_put_vec(yvals);
            //     cout<<"r = "<<r<<endl;
            //     cout<<"dr = "<<dr<<endl;
    
            step(rhs, y, r_0, dr0 ); 

            if(r_0<rf)
            {
                yvals = y;
                r = r_0;
                dr = dr0;
                check_field_isfinite(r,yvals,"solution is nan");

                

            }
            else{
                counter = 1;
            }
            
        }
            dr = rf - r;
            step(rhs, yvals, r, dr );
            check_field_isfinite(r,yvals,"solution is nan");

            // out_put_vec(yvals);
            //     cout<<"r = "<<r<<endl;
            //     cout<<"dr = "<<dr<<endl;
        
    }
    else if(side == "right"){
        assert(r0<R_star);
        assert(r0>rf);
        if(dr>0){ dr = -dr; }
        if(get_series)
        {   
            
            
            assert(ics.size()==3);
            series_surf( r0, ics, series_sol) ;

            // cout<<"get_series_true"<<endl;
            // cout<<side<<endl;

        }
            
        
        yvals = series_sol;
        // cout<<"part = "<<sol_type<<endl;
        // out_put_vec(yvals);
        // cout<<"series at surface"<<endl;
        // out_put_vec(yvals);
        // cout<<"---"<<endl;
        double r = r0;
        int counter = 0;
        while(counter ==0)
        {   
            vector<double> y = yvals;
            double r_0 = r;
            double dr0 = dr;
            
            step(rhs, y, r_0, dr0 );
            
            if(r_0>rf)
            {
                yvals = y;
                r = r_0;
                dr = dr0;
                check_field_isfinite(r,yvals,"solution is nan");

                // out_put_vec(yvals);
                // cout<<"r = "<<r<<endl;
                // cout<<"dr = "<<dr<<endl;

            }
            else{
                counter = 1;
            }
            
        }
            dr = rf - r;
            step(rhs, yvals, r, dr );
            check_field_isfinite(r,yvals,"solution is nan");

            // out_put_vec(yvals);
            //     cout<<"r = "<<r<<endl;
            //     cout<<"dr = "<<dr<<endl;

    }
    // yfinal = yvals;
}