#ifndef _OUTPUT_FILES_HPP_
#define _OUTPUT_FILES_HPP_

#include <string>
#include <vector>

 void write_vec(std::string path, const std::vector<double> vecl, const std::string name_v);
 void write_double(std::string path, double omega, const std::string name_v);
 void check_field_isfinite(
    const double xi,
    const std::vector<double> vec, std::string message);

void out_put_vec(std::vector<double> vec);

void combine_vec(const std::vector<double> vecl, const std::vector<double> vecr, std::vector<double> &vec);

double fzeta_func(double theta, double n);


double feta_func(double theta, double n);

double feta_der_func(double theta, double n);

double feta_der_2_func(double theta, double n);

void subdivide(double start, double stop, int N, std::vector<double> &vec);

void write_vec_omega_k2_tau(std::string path, 
const std::vector<double> Omega_vals, 
const std::vector<double> K2, 
const std::vector<double> K2tau, const std::string name_v);

void write_vec_omega_k2(std::string path, 
const std::vector<double> Omega_vals, 
const std::vector<double> K2,  const std::string name_v);

void shift_vec(double val, std::vector<double> &vec);


void trapezoid_integrate(
  double n,
  std::string type_visc,
  std::vector<double> thetavals, 
  std::vector<double> lambdavals,
  std::vector<double> xi, double &dens_avg, double &visc_avg);

// double get_random();


#endif