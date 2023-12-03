#include <cassert>
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <iomanip>
#include <fstream>
#include<cassert>

#include "sim_params.hpp"
/*===========================================================================*/
string Sim_params::read_sim_params(const string output_dir, const string find_this_var)
{
   ifstream infile(output_dir+"/sim_params.txt");
   assert(infile.good());

   string name, val;

   while (infile >> name >> val) {
      if (name==find_this_var) {
         return val;
      }
   }
   cout << "ERROR: did not find " << find_this_var << endl;
   std::exit(0);
}
/*===========================================================================*/
Sim_params::Sim_params(const string output_dir)
{  
   n= stod(read_sim_params(output_dir,"n"));
   N1= stod(read_sim_params(output_dir,"N1"));
   b= stod(read_sim_params(output_dir,"b"));
   xil= stod(read_sim_params(output_dir,"xil"));
   N_interp_LE= stoi(read_sim_params(output_dir,"N_interp_LE"));
   steps_num = stoi(read_sim_params(output_dir,"steps_num"));
   N_write = stoi(read_sim_params(output_dir,"N_write"));
//----------------------------------------------------------------------
   Omegalow = stod(read_sim_params(output_dir,"Omegalow"));
   Omegahigh = stod(read_sim_params(output_dir,"Omegahigh"));
   div_Omega = stoi(read_sim_params(output_dir,"div_Omega"));
   
   factor_match = stod(read_sim_params(output_dir,"factor_match"));

   // solve_visc = stoi(read_sim_params(output_dir,"solve_visc"));
   // visc_type = read_sim_params(output_dir,"visc_type");

}
/*============================================================================*/
Sim_params::~Sim_params(void)
{
}
/*============================================================================*/
void Sim_params::output_parameters()
{
//============================================================================
   cout<<"n = "<<n<<endl;
   cout<<"b = "<<b<<endl;
   cout<<"xil = "<<xil<<endl;
   // cout<<"xir = "<<xir<<endl;
//============================================================================
//============================================================================
}

/*============================================================================*/
