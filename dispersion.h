#ifndef DISPERSION_H_
#define DISPERSION_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

int whatsthatsign(double x);

struct dispParams{ double T;double M;double k0; double k1;};

double LpB(double p, double M, double k0, double k1);
double LmB(double p, double M, double k0, double k1);

double TRKS_integrand(double p, void *d);
double TRKS(double T, double M,double k0, double k1);
double TRUS_integrand(double p, void *d);
double TRUS(double T, double M,double k0, double k1);

double sigmaA(struct dispParams pp);
double sigmaB(struct dispParams pp);
double dispEquation(double w1, void * d);
double dispSolved(double w0, struct dispParams p);

#endif

