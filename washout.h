#ifndef WASHOUT_H_
#define WASHOUT_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

struct washoutParams{ double T;double M;};
double washout_integrand (double x,  void * p);
double washout(double eta,double T, double M);

#endif

