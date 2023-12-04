#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <fstream>
#include <cmath>
using std::pow;
using std::abs;
using std::log;
using std::isfinite;
#include<cassert>
// #include <random>

#include "outputfiles.hpp"

void write_vec(std::string path, const std::vector<double> vecl, const std::vector<double> vecr, const string name_v){
  string name = path + name_v + ".dat";
  std::ofstream write_output(name , std::ios::app);
  assert(write_output.is_open());
  write_output.precision(16);

  int nx = vecl.size();
  for(int i = 0; i<nx; ++i ){
    write_output<<vecl[i]<<",";
  }
  for(int i = 0; i<nx-2; ++i ){
    write_output<<vecr[i]<<",";
  }
  write_output<<vecr[nx-2]<<"\n";

  write_output.flush();
  write_output.close();

}
void write_double(std::string path,double omega, const string name_v){
  string name = path + name_v + ".dat";
  std::ofstream write_output(name , std::ios::app);
  assert(write_output.is_open());
  write_output.precision(16);

  write_output<<omega<<",";

  write_output.flush();
  write_output.close();

}

void check_field_isfinite(
    const double xi,
    const std::vector<double> vec, std::string message)
{

    size_t index=0;
    for (auto x: vec) {
       if (!isfinite(x)) {
          cout<<"NaN;  ";
          cout<<"index =  "<<index<<" xi = "<<xi<<endl;
          cout<<message<<endl;
          std::exit(0);
       }
       index++;
    }
}

void out_put_vec(std::vector<double> vec)
{
  int len = vec.size();
  for(int i=0;i<len; ++i)
  {
    cout<<"i = "<<i<<"; val = "<<vec[i]<<endl;
  }
}

void combine_vec(const std::vector<double> vecl, const std::vector<double> vecr, std::vector<double> &vec)
{
  int N = vecl.size();

  for(int i=0; i<N; ++i)
  {
    vec[i] = vecl[i];
    // cout<<i<<endl;
    // cout<<vec[i]<<endl;
  }
  for(int i=0; i<N-1; ++i)
  {
    vec[N+i] = vecr[i];
    // cout<<N+i<<endl;
    // cout<<vec[i]<<endl;
  }
}

double fzeta_func(double theta, double n)
{
  return pow(theta,n+2);
}




double feta_func(double theta, double n)
{
 double i = 4;
  return pow(theta,n+i);
  // return 0.;

}

double feta_der_func(double theta, double n)
{

double i = 4;
return (n+i)*pow(theta,n+i-1);
  // return 0.;
}

double feta_der_2_func(double theta, double n)
{
double i = 4;
return (n+i)*(n+i-1)*pow(theta,n+i-2);
// return 0.;
}







void subdivide(double start, double stop, int N, std::vector<double> &vec)
{
  {
    int Nsize = vec.size();
    assert(Nsize==N);
  }
  
  double dx = (stop - start)/(N-1);

  for(int i=0; i<N; ++i)
  {
    vec[i] = start + i*dx;
  }
}

void write_vec_omega_k2_tau(std::string path, 
const std::vector<double> Omega_vals, 
const std::vector<double> K2, 
const std::vector<double> K2tau, const std::string name_v)
{
   string name = path+"/" + name_v + ".dat";
    std::ofstream write_output(name , std::ios::app);
    assert(write_output.is_open());
  write_output.precision(16);

  int nx = Omega_vals.size();
  for(int i = 0; i<nx; ++i ){
    write_output<<Omega_vals[i]<<"\t"<<K2[i]<<"\t"<<K2tau[i]<<"\n";
  }
 

  write_output.flush();
  write_output.close();
}

void write_vec_omega_k2(std::string path, 
const std::vector<double> Omega_vals, 
const std::vector<double> K2,  const std::string name_v)
{
  string name = path+"/" + name_v + ".dat";
    std::ofstream write_output(name , std::ios::app);
    assert(write_output.is_open());
  write_output.precision(16);

  int nx = Omega_vals.size();
  for(int i = 0; i<nx; ++i ){
    write_output<<Omega_vals[i]<<"\t"<<K2[i]<<"\n";
  }
 

  write_output.flush();
  write_output.close();
}

void shift_vec(double val,std::vector<double> &vec)
{
  int N = vec.size();

  for(int i=0; i<N; ++i)
  {
    vec[i] -= val;
  }
}

void trapezoid_integrate(
  double n,
  std::string type_visc,
  std::vector<double> thetavals, 
  std::vector<double> lambdavals,
  std::vector<double> xi, double &dens_avg, double &visc_avg)
{
  dens_avg = 0.; 
  visc_avg = 0.;
  double avg = 0.;
  int N  = xi.size();
  for(int i=1; i<N; ++i)
  { 
    double thetaim1 =thetavals[i-1] ;
    double thetai = thetavals[i];
    double xii = xi[i];
    double xiim1 = xi[i-1];
    double lambdai = lambdavals[i];
    double lambdaim1 = lambdavals[i-1];
    double fim1 = pow(thetaim1, n)*pow(xiim1,2.)*exp(lambdaim1/2.);
    double fi = pow(thetai, n)*pow(xii,2.)*exp(lambdai/2.);
    double dx = xii-xiim1;

    dens_avg += 0.5*(dx)*(fim1 + fi);

    avg += 0.5*dx*(pow(xiim1,2.)*exp(lambdaim1/2.) + pow(xii,2.)*exp(lambdai/2.) );

    if(type_visc == "bulk")
    {
      fim1 = fzeta_func(thetaim1,n)*pow(xiim1,2.)*exp(lambdaim1/2.);
      fi = fzeta_func(thetai, n)*pow(xii,2.)*exp(lambdai/2.);

    }
    else{
     fim1 = feta_func(thetaim1,n)*pow(xiim1,2.)*exp(lambdaim1/2.);
      fi = feta_func(thetai, n)*pow(xii,2.)*exp(lambdai/2.);
    }
    visc_avg += 0.5*(dx)*(fim1 + fi);
  }
  dens_avg = dens_avg/avg;
  visc_avg = visc_avg/avg;

}