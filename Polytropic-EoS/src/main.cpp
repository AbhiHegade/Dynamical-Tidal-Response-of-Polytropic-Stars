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
<<xi1_true<<"\t"
<<Ca<<endl;
write_bkg.close();

// /*-----------------------------------------------------------*/
// /*-----------------------------------------------------------*/
int divisions = sp.div_Omega;
assert(sp.Omegahigh>sp.Omegalow);
double deltaomega = (sp.Omegahigh - sp.Omegalow)/(divisions-1);
// vector<double> ks(divisions);
// vector<double> omegas(divisions);
for(int i=0; i<divisions; ++i)
{  
  cout<<"Start Run"<<endl;
  cout<<"===================="<<endl;
  double omega = sp.Omegalow +i*deltaomega;
  double k2 = 0.;
//   double k2tau = 0.;
  cout<<"i = "<<i<<endl;
  cout<<"omega = "<<omega<<endl;
  double l=2.;
  double gamma_0 = 1.0 + 1.0/sp.N1;
    
  CONSERVATIVE_TIDE pf(l, omega, sp.n, sp.b, le.nu0, gamma_0, muf,xi1_true, theta, lambda, nu);
  
  double dr_internal = 1e-10;

  int N_write = sp.N_write;

  vector<double> icsl(4), icsr(4);

  vector<double> H0_v(N_write);
  vector<double> W0_v(N_write);
  vector<double> V_v(N_write);
  vector<double> H1_v(N_write);
  vector<double> r_v(N_write);

  vector<double> vals_origin(2);

  double normk2 =0.; 
    
  pf.solve_and_iterate_and_write(N_write, sp.steps_num, dr_internal, 
    sp.xil, sp.factor_match*(sp.xil + xi1_true), xi1guess, 
    k2, H0_v, W0_v, V_v, H1_v, r_v, vals_origin, normk2);



    write_output<<omega<<"\t"<<k2<<"\t"
        <<normk2<<"\n";
        write_output.flush();
        cout<<"End Run."<<endl;
        cout<<"===================="<<endl;

}

}

//=======================================================
//=======================================================
void get_diff_radius(double n, double brel, double thetaf, double muR0, double xi1, double &xitrue )
{ 
    
    double a = -(muR0/((2*brel*muR0*(1 + n) - xi1)*xi1));

    double diff = thetaf/(a);

    xitrue = xi1 + diff;
 
    // cout<<"diff = "<<diff<<endl;
}