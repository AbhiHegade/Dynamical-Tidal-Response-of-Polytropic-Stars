#include<iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <sstream>
#include <string>
using std::string;
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;
#include <boost/math/interpolators/barycentric_rational.hpp>
using boost::math::interpolators::barycentric_rational;
#include <cmath>
using std::pow;
using std::abs;
using std::log;
#include <random>
#include <fstream>
#include <cassert>
//---------------------------------------
#include "lane-emden-solver.hpp"
#include "sim_params.hpp"
#include "outputfiles.hpp"
#include "conservative_tides.hpp"
#include "dissipative_tides_bulk.hpp"
#include "shear_funcs.hpp"

void get_diff_radius(double n ,double brel, double thetaf, double muR0, double xi1, double &xitrue );
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
int main(int argc, char *argv[]) {

assert(argc == 2);
string path = argv[1];
Sim_params sp(path);

string name = path+"/" + "vals" + ".dat";
std::ofstream write_output(name , std::ios::app);
assert(write_output.is_open());
write_output.precision(16);
//----------------------------------------------------------------------
//--------------Solving Relativistic Lane Emden Equation Now------------
LE_SOLVER le(sp.n, sp.b);

double xi1guess = 0.,thetaf = 0., muf = 0., nuf=0.;
int N = 100000;
le.solve_adaptive(sp.xil,N,xi1guess, thetaf, muf, nuf);

double xi1_true =  xi1guess;

get_diff_radius(sp.n, sp.b, thetaf, muf, xi1guess, xi1_true);

cout<<"xi1_true = "<<xi1_true<<endl;

N = sp.N_interp_LE;
vector<double> xivals_LE_v(N), theta_v(N), lambda_v(N), nu_v(N);

double dxi_theta_mu = 0.;

le.solve_and_interp(N, sp.xil, xi1guess, dxi_theta_mu, thetaf, muf, nuf, xivals_LE_v, theta_v, lambda_v , nu_v);


cardinal_cubic_b_spline<double> theta(theta_v.begin(), theta_v.end(), sp.xil, dxi_theta_mu);
cardinal_cubic_b_spline<double> lambda(lambda_v.begin(), lambda_v.end(), sp.xil, dxi_theta_mu); 
cardinal_cubic_b_spline<double> nu(nu_v.begin(), nu_v.end(), sp.xil, dxi_theta_mu); 

double Ca = sp.b*(sp.n+1)*muf/xi1_true;
std::ofstream write_bkg(path+"/" + "bkg" + ".dat", std::ios::app);
assert(write_bkg.is_open());
write_bkg.precision(16);
write_bkg<<sp.n<<"\t"
<<sp.b<<"\t"
<<muf<<"\t"
<<muf*pow(sp.b,(3.-sp.n)/2.)<<"\t"
<<xi1_true<<"\t"
<<xi1_true*pow(sp.b, (1-sp.n)/2.)<<"\t"
<<Ca<<endl;
write_bkg.close();

if(sp.solve_tides==0)
{
 std::exit(0);
}

// /*-----------------------------------------------------------*/
double dens_avg = 0.;
double visc_bulk_avg = 0.;
double visc_shear_avg=0.;
if(sp.solve_visc)
{
    // Compute density and viscosity average

    trapezoid_integrate(sp.n, 
    "bulk", theta_v, lambda_v,  xivals_LE_v, dens_avg, visc_bulk_avg);

    trapezoid_integrate(sp.n, 
    "shear", theta_v, lambda_v,  xivals_LE_v, dens_avg, visc_shear_avg);

    cout<<"dens_avg = "<<dens_avg<<endl;
    cout<<"visc_bulk_avg = "<<visc_bulk_avg<<endl;
    cout<<"visc_shear_avg = "<<visc_shear_avg<<endl;
    cout<<"dens_avg/visc_bulk_avg = "<<dens_avg/visc_bulk_avg<<endl;
}
// /*-----------------------------------------------------------*/
int divisions = sp.div_Omega;
assert(sp.Omegahigh>sp.Omegalow);
double deltaomega = (sp.Omegahigh - sp.Omegalow)/(divisions);
for(int i=0; i<divisions+1; ++i)
{  
  cout<<"Start Run"<<endl;
  cout<<"===================="<<endl;
  double omega = sp.Omegalow +i*deltaomega;
  double k2 = 0.;
  double k2tau = 0.;
  cout<<"i = "<<i<<endl;
  cout<<"omega = "<<omega<<endl;
  double l=2.;
  double gamma_0 = 1.0 + 1.0/sp.N1;
    
  CONSERVATIVE_TIDE pf(l, omega, sp.n, sp.b, le.nu0, gamma_0, muf,xi1_true, theta, lambda, nu);
  
  double dr_internal = 1e-10;

  int N_write = sp.N_write;

  vector<double> H0_v(N_write);
  vector<double> W0_v(N_write);
  vector<double> V_v(N_write);
  vector<double> H1_v(N_write);
  vector<double> xi_v(N_write);

  vector<double> vals_origin(2);

  double normk2 =0.; 

  if(sp.solve_visc || sp.write_vec)
  {
    

    
    pf.solve_and_iterate_and_write(N_write, sp.steps_num, dr_internal, 
      sp.xil, sp.factor_match*(sp.xil + xi1_true), xi1guess, 
      k2, H0_v, W0_v, V_v, H1_v, xi_v, vals_origin, normk2);
    
    if(sp.write_vec)
    {
      write_vec(path + "/", H0_v, "H0" );
      write_vec(path + "/", W0_v, "W0");
      write_vec(path + "/", V_v, "V");
      write_vec(path + "/", H1_v, "H1");
      if(i==0)
      {
        write_vec(path + "/", xi_v,"xi");
      }

    }
  }
  else{
    vector<double> icsl(4), icsr(4);
    pf.solve_and_iterate(sp.steps_num, dr_internal, 
      sp.xil, sp.factor_match*(sp.xil + xi1_true), xi1guess, 
      k2, normk2, icsl, icsr);
  }
  

  
  if(sp.solve_visc)
  {
    cardinal_cubic_b_spline<double> H0_func(H0_v.begin(), H0_v.end(), sp.xil, xi_v[1]-xi_v[0]);
    cardinal_cubic_b_spline<double> W0_func(W0_v.begin(), W0_v.end(), sp.xil, xi_v[1]- xi_v[0]);
    cardinal_cubic_b_spline<double> V_func(V_v.begin(), V_v.end(), sp.xil, xi_v[1]- xi_v[0]);
    cardinal_cubic_b_spline<double> H0der_func(H1_v.begin(), H1_v.end(), sp.xil, xi_v[1]- xi_v[0]);

    cout<<"=============================================="<<endl;
    cout<<"Solving perturbed problem now"<<endl;


    DISS_BULK_TIDE pfbulk(l, omega, sp.n, sp.b, le.nu0, gamma_0, muf, xi1_true,
    theta, lambda, nu,
    fzeta_func,
    H0_func, 
    W0_func, 
    V_func,
    vals_origin[0], vals_origin[1]);


    dr_internal = 1e-8;

    double normk2tau = 0.;

    pfbulk.solve_and_iterate(sp.steps_num, dr_internal, sp.xil, sp.factor_match*(sp.xil + xi1_true),xi1guess, k2tau, normk2tau );

    double p2 = (muf/xi1_true)*(k2tau/k2)*(sp.n+1)*(dens_avg/visc_bulk_avg);

    double hat_tau = k2tau/k2;

    write_output
    <<omega<<"\t"<<k2<<"\t"<<normk2<<"\t"
    <<hat_tau<<"\t"<<p2<<"\t"<<normk2tau<<"\t";
  
//------------------------------------------------------------------------
    cout<<"Solving for shear lag"<<endl;

    DISS_SHEAR_TIDE pfshear(l, omega, sp.n, sp.b, le.nu0, gamma_0, muf, xi1_true, theta, lambda, nu,
    feta_func, feta_der_func, feta_der_2_func, 
    H0_func,W0_func, V_func, H0der_func,
    vals_origin[0], vals_origin[1]);

    dr_internal = 1e-8;
    normk2tau = 0.;
    k2tau = 0.;

    out_put_vec(vals_origin);

    // cout<<"here"<<endl;

    pfshear.solve_and_iterate(sp.steps_num, dr_internal, sp.xil, sp.factor_match*(sp.xil + xi1_true),xi1guess, k2tau, normk2tau );

    // cout<<"here"<<endl;

    p2 = (muf/xi1_true)*(k2tau/k2)*(sp.n+1)*(dens_avg/visc_shear_avg);

    hat_tau = k2tau/k2;
    
    write_output<<hat_tau<<"\t"<<p2<<"\t"<<normk2tau<<"\n";
    write_output.flush();
    cout<<"End Run."<<endl;
    cout<<"===================="<<endl;

  }
  else{


    write_output<<omega<<"\t"<<k2<<"\t"
    <<normk2<<"\n";
    write_output.flush();
    cout<<"End Run."<<endl;
    cout<<"===================="<<endl;
  }

}

}

//=======================================================
//=======================================================
void get_diff_radius(double n, double brel, double thetaf, double muR0, double xi1, double &xitrue )
{ 
    
    double a = -(muR0/((2*brel*muR0*(1 + n) - xi1)*xi1));

    double diff = thetaf/(a);

    xitrue = xi1 + diff;
 
}
