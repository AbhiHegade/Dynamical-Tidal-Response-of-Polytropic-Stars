#ifndef _READ_PARAM_FILE_HPP_
#define _READ_PARAM_FILE_HPP_

#include<string>

class Sim_params
{
public:
  double n;
  double N1;
  double b;
  double xil;
  int N_interp_LE;
  int steps_num;
  int N_write;
  int solve_tides;
  //-----------
  double Omegalow;
  double Omegahigh;
  int div_Omega;
  double factor_match;
  int solve_visc;
  int write_vec;
  // std::string visc_type;
  //-----------
//------------------------------------------------------------------------------
  Sim_params(const std::string output_dir);
//------------------------------------------------------------------------------
  ~Sim_params(void);
//------------------------------------------------------------------------------
  void output_parameters();
private:
   std::string read_sim_params(const std::string output_dir, const std::string find_this_var);

};















#endif
