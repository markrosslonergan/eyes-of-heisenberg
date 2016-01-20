#ifndef DISPERSION_H_
#define DISPERSION_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

int whatsthatsign(double x);


struct dispParams{ double T;double M;double k0;double k1;bool use_nf;};
struct dispParams2{ double T;double M;double k1;bool use_c; bool use_nf;};

double nB(double p, double T);
double nF(double p, double T, double M);

double LpB_negative(double p, double M, double k0, double k1);

double LpB(double p, double M, double k0, double k1);
double LmB(double p, double M, double k0, double k1);
double LpF(double p, double M, double k0, double k1);
double LmF(double p, double M, double k0, double k1);
double LpB_pole_top(double M, double k0, double k1);
double LpB_pole_bot(double M, double k0, double k1);



double TRKS_integrand(double p, void *d);
double TRKS(double T, double M,double k0, double k1,bool use_nf);
double TRUS_integrand(double p, void *d);
double TRUS(double T, double M,double k0, double k1,bool use_nf);
double TRS_integrand(double p, void *d);
double TRS(double T,double M,double k0,double k1,bool use_nf);

double my_log_test();

double sigmaA(double T, double M, double k0, double k1, bool use_nf);
double sigmaB(double T, double M, double k0, double k1, bool use_nf);
double sigmaC(double T, double M, double k0, double k1, bool use_nf);
double dispEquation(double k0, void * d);
double dispEquationN(double k0, void * d);
double dispSolved(double T, double M, double k1,bool use_c, bool use_nf);

double plot_disp(double M, double T, bool use_c, bool use_nf, double k1oT_min, double k1oT_max);


#endif

