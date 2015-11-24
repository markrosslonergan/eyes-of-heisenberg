#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <nlopt.hpp>
#include <iomanip>
#include <sys/time.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_log.h>
//#include <gsl/pow.h>

#include "dispersion.h"

#define GW 2.0
#define PI 3.14149
#define REL 0.01
#define ABS 0

double washout_integrand (double x,  void * p) {
	struct washoutParams * params  = (struct washoutParams *)p;
	double M = params->M;
	double T = params->T;
	double a=1.0;

	double top =exp(a*a*M*M/(4*T*T*x)+x)+1;
	double bottom =exp(a*a*M*M/(4*T*T*x)+x)-exp(x);

	return  exp(x)/((exp(x)+1.0)*(exp(x)+1.0))*gsl_sf_log(top/bottom);
}


double washout(double eta,double T, double M)
{
	double Y1=1.0;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	struct washoutParams params = {T,M};
	double result,error;
	gsl_function F;
	F.function = &washout_integrand;
	F.params = &params;
	
	gsl_integration_qags(&F,0.0, 50, ABS, REL, 5000, w, &result, &error); 
	gsl_integration_workspace_free (w);

	return -3*Y1*Y1*M*M*eta/(8*PI*PI*PI*T)*result;

}


