#ifndef EPSILON_H_
#define EPSILON_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

struct epsilonParams { double T; double M; };

double I0(double y0, double y);
double I1(double y0, double y);

double sigma0N(double T, double k0, double k);
double sigmaiN(double T, double M, double k0, double k);

double delta_FN1(double T, double M, double k);

double epsilon_numerator_integrand (double k,  void * p); 
double epsilon_denominator_integrand (double k,  void * p); 
double epsilon_numerator(double T, double M, double cutoff);
double epsilon_denominator(double T, double M, double cutoff);

double epsilon_eff_ratio (double T, double M, double cutoff);



#endif

