#ifndef _LOCAL_EXPANSIONS_HPP_
#define _LOCAL_EXPANSIONS_HPP_

#include <vector>
#include <Eigen/Dense>

void get_matrix_ord_2_origin(double l, double Omega, double nu0, double n, double brel, 
double gamma_0,
double xi, Eigen::MatrixXd &A);

void get_matrix_ord_2_origin_just_H_W(double l, double Omega, double nu0, double n, double brel, 
double gamma_0,
double xi, Eigen::MatrixXd &A);
//----------------------------------------------
//----------------------------------------------
class SERIES_ORIGIN_SOURCELESS{

    public:
    SERIES_ORIGIN_SOURCELESS(double l, double Omega, double nu0, double n, double brel, double gamma_0);
    ~SERIES_ORIGIN_SOURCELESS(void);
//----------------------------------------------
    double l; 
    double Omega; 
    double nu0;
    double n; 
    double brel; 
    double gamma_0;

    void operator()( double xi, const std::vector<double> yin, std::vector<double> &y );


    private:


};

//----------------------------------------------
class SERIES_SURF{

    public:
    SERIES_SURF(double l, double Omega, double n, double brel, double gamma_0, double muR0, double xi1);
    ~SERIES_SURF(void);
//----------------------------------------------
   double l; 
   double Omega;
   double n; 
   double brel; 
   double gamma_0; 
   double muR0; 
   double xi1;

    void operator()( double xi, const std::vector<double> yin, std::vector<double> &y );


    private:


};
//------------------------------------------------------------------
void obtain_H0_and_H1_ext(double l, double Omega, double n, double brel, double mu1, double xi1, std::vector<double> &H0_ext, std::vector<double> &H1_ext );
//------------------------------------------------------------------














#endif