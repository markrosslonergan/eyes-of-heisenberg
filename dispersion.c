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
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
//#include <gsl/pow.h>

#include "dispersion.h"

#define GW 2.0
#define PI 3.14149
#define REL 0.0001
#define ABS 0

double nB(double p, double T){

return 1.0/(exp(p/T)-1.0);
}


int whatsthatsign(double x){
	if (x > 0) return 1;
	if (x < 0) return -1;
return 0;
}

double LpB(double p, double M, double k0, double k1){

	double top = pow(-k1*k1+k0*k0-M*M,2)+4*k1*(-k1*k1+k0*k0-M*M)*p - 4*(-k1*k1+k0*k0)*p*p;
	double bot = pow(-k1*k1+k0*k0-M*M,2)-4*k1*(-k1*k1+k0*k0-M*M)*p - 4*(-k1*k1+k0*k0)*p*p;
	return gsl_sf_log_abs(top/bot);
}

double LmB(double p, double M, double k0, double k1){
	double top = pow(-k0*k0+k1*k1+M*M,2)-4*pow(k0+k1,2)*p*p;
	double bot = pow(-k0*k0+k1*k1+M*M,2)-4*pow(k0-k1,2)*p*p;

	return gsl_sf_log_abs(top/bot);
}



double TRKS_integrand(double p, void *d){
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k1 = params->k1;
	double k0 = params->k0;
	double T = params->T;

	return nB(p,T)/(2*2*PI*2*PI)*p*(4+(M*M-k1*k1+k0*k0)*LpB(p,M,k0,k1)/(2*p*k1));
}


double TRKS(double T, double M,double k0, double k1)
{

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	struct dispParams params = {T,M,k0,k1};
	double result,error;
	gsl_function F;
	F.function = &TRKS_integrand;
	F.params = &params;
	
	gsl_integration_qags(&F,0.0, 30*T, ABS, REL, 5000, w, &result, &error); 
	gsl_integration_workspace_free (w);

	return result;

}


double TRUS_integrand(double p, void *d){
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k1 = params->k1;
	double k0 = params->k0;
	double T = params->T;

	return nB(p,T)/(2*2*PI*2*PI*k1)*(p*LmB(p,M,k0,k1)+k0*LpB(p,M,k0,k1));
}


double TRUS(double T, double M,double k0, double k1)
{

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	struct dispParams params = {T,M,k0,k1};
	double result,error;
	gsl_function F;
	F.function = &TRUS_integrand;
	F.params = &params;
	
	gsl_integration_qags(&F,0.0, 30*T, ABS, REL, 5000, w, &result, &error); 
	gsl_integration_workspace_free (w);

	return result;

}
double sigmaA(struct dispParams p){
	return 1/(p.k1*p.k1)*(TRKS(p.T,p.M,p.k0,p.k1)-p.k0*TRUS(p.T,p.M,p.k0,p.k1));
}
double sigmaB(struct dispParams p){
	return (pow(p.k0/p.k1,2)-1)*TRUS(p.T,p.M,p.k0,p.k1)-(p.k0/pow(p.k1,2))*TRKS(p.T,p.M,p.k0,p.k1);
}

double dispEquation(double w1, void * d){
	struct dispParams * params = (struct dispParams *)d;
	double T = params->T;
	double M = params->M;
	double k0 = params->k0;
	params->k1=w1;

	return (1-sigmaA(*params))*k0-sigmaB(*params)-(1-sigmaA(*params))*w1;

}

double dispSolved(double w0, struct dispParams p){

	int iter = 0, max_iter =200;
	int status;
	double r = 0;

	const gsl_root_fsolver_type * T   = gsl_root_fsolver_brent;
	gsl_root_fsolver * s     = gsl_root_fsolver_alloc (T);
	//const gsl_root_fdfsolver_type * T   = gsl_root_fdfsolver_newton;
	//gsl_root_fdfsolver * s     = gsl_root_fdfsolver_alloc (T);


	p.k0=w0;

	
	double x_lo=w0;
	double x_hi=w0;
	
	double multiplier = 0.99;
	double m = 1;
	for(m = 1; whatsthatsign(dispEquation(w0*(m*multiplier),&p))==whatsthatsign(dispEquation(w0/(m*multiplier),&p)); m=m+1)
	{
		//std::cout<<m*multiplier<<" "<<dispEquation(w0*(m*multiplier),&p)<<" "<<dispEquation(w0/(m*multiplier),&p)<<std::endl;
	
	}
		//std::cout<<m*multiplier<<" "<<dispEquation(w0*(m*multiplier),&p)<<" "<<dispEquation(w0/(m*multiplier),&p)<<std::endl;
		x_lo=std::min(w0*m*multiplier,w0/(m*multiplier));
		x_hi=std::max(w0*m*multiplier,w0/(m*multiplier));

	gsl_function F;
	F.function =&dispEquation;
	F.params =&p;

	gsl_root_fsolver_set(s,&F,x_lo,x_hi);
	


	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_x_lower(s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo,x_hi,0,0.0001);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);

	return r;

}

