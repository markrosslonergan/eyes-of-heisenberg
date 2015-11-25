#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <sys/time.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_log.h>
//#include <gsl/pow.h>

#include "epsilon.h"

#define GW 2.0
#define PI 3.14149
#define REL 0.01
#define ABS 0

double I0(double y0, double y)
{
	return gsl_sf_log((exp(y0+y)-1.0)/(exp(y0-y)-1.0))-y;
}

double I1(double y0, double y)
{
	double p1 = 0.5*(y0+y)*gsl_sf_log((1.0+exp((y0+y)/2.0))/(1-exp((-y0+y)/2))); 
	double p2 = -0.5*(y0-y)*gsl_sf_log((1.0+exp((y0-y)/2.0))/(1-exp((-y0-y)/2))); 
	return p1+p2+gsl_sf_dilog(-exp(y0/2.0+y/2.0))+gsl_sf_dilog(exp(-y0/2.0-y/2.0))-gsl_sf_dilog(-exp(y0/2.0-y/2.0))-gsl_sf_dilog(exp(-y0/2.0+y/2.0));
}

double sigma0N(double T, double k0, double k)
{
	return (GW*T*T)/(8*PI*k)*I1(k0/T,k/T);
}

double sigmaiN(double T, double M, double k0, double k)
{
	return GW*T*T/(8*PI*k)*(k0/k*I1(k0/T,k/T)-M*M/(2*k*T)*I0(k0/T,k/T))*k/k;
}

double delta_FN1(double T, double M, double k)
{
	double beta = 1.0/T;

	return 1.0/(exp(beta*sqrt(k*k+M*M))+1);
}


double epsilon_numerator_integrand (double k,  void * p) {
	struct epsilonParams * params  = (struct epsilonParams *)p;
	double M = params->M;
	double T = params->T;
	double k0 = sqrt(M*M+k*k);


	return k*k/k0*(pow(sigma0N(T,k0,k),2.0)-pow(sigmaiN(T,M,k0,k),2.0))*delta_FN1(T,M,k); 
}

double epsilon_denominator_integrand (double k,  void * p) {
	struct epsilonParams * params  = (struct epsilonParams *)p;
	double M = params->M;
	double T = params->T;
	double k0 = sqrt(M*M+k*k);


	return k*k/k0*(k0*sigma0N(T,k0,k)-k*sigmaiN(T,M,k0,k))*delta_FN1(T,M,k); 
}

double epsilon_numerator(double T, double M, double cutoff)
{
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	struct epsilonParams params = {T,M};
	double result,error;
	gsl_function F;
	F.function = &epsilon_numerator_integrand;
	F.params = &params;
	
	gsl_integration_qags(&F,cutoff ,100*T,ABS, REL, 5000, w, &result, &error); 
	gsl_integration_workspace_free (w);

	return result;

}

double epsilon_denominator(double T, double M, double cutoff)
{
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	struct epsilonParams params = {T,M};
	double result,error;
	gsl_function F;
	F.function = &epsilon_denominator_integrand;
	F.params = &params;
	
	gsl_integration_qags(&F,cutoff,100*T,ABS, REL, 5000, w, &result, &error); 
	gsl_integration_workspace_free (w);

	return result;

}

double epsilon_eff_ratio (double T, double M, double cutoff)
{
return 16.0*PI/GW*epsilon_numerator(T,M,cutoff)/epsilon_denominator(T,M,cutoff);
}
