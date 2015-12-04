#ifndef DISPERSION_H_
#define DISPERSION_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

int whatsthatsign(double x);


struct dispParams{ double T;double M;double k1;double k0;};
struct dispParams2{ double T;double M;double k1;};

double nB(double p, double T);
double nF(double p, double T, double M);


double LpB(double p, double M, double k0, double k1);
double LmB(double p, double M, double k0, double k1);
double LpF(double p, double M, double k0, double k1);
double LmF(double p, double M, double k0, double k1);


double TRKS_integrand(double p, void *d);
double TRKS(double T, double M,double k0, double k1);
double TRUS_integrand(double p, void *d);
double TRUS(double T, double M,double k0, double k1);
double TRS_integrand(double p, void *d);
double TRS(double T,double M,double k0,double k1);


double sigmaA(double T, double M, double k0, double k1);
double sigmaB(double T, double M, double k0, double k1);
double sigmaC(double T, double M, double k0, double k1);
double dispEquation(double k0, void * d);
double dispSolved(double T, double M, double k1);

std::vector<double > plotDispersion(double M,double T,double kmin, double kmax);

#endif

