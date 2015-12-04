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
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "dispersion.h"

#define GW 2.0
#define PI 3.14149
//#define REL 0.000001//0.0001
#define REL 1.0e-6
#define ABS 0.0
//#define ABS 0
#define DEBUG_MODE false
//#define T_TRESHOLD 1e-10
#define T_TRESHOLD 1e-10



double nB(double p, double T){
   double ans = 0.0;
    if (p < 1e8*T){ ans = 1.0/(exp(p/T)-1.0);}
    
    return ans;
}



double nF(double p, double T, double M){
    double ans = 0.0;
    //if (sqrt(p*p+M*M) < 1e5*T){ ans = 1.0/(exp(sqrt(p*p+M*M)/ T)+1.0);}
    ans =1.0/(exp(sqrt(p*p+M*M)/T)+1.0); 
    return ans;

}

int whatsthatsign(double x){
	if (x > 0) return 1;
	if (x < 0) return -1;
return 0;
}




double LpB(double p, double M, double k0, double k1){

	double top = pow(-k1*k1+k0*k0-M*M,2)+4*k1*(-k1*k1+k0*k0-M*M)*p - 4*(-k1*k1+k0*k0)*p*p;
	double bot = pow(-k1*k1+k0*k0-M*M,2)-4*k1*(-k1*k1+k0*k0-M*M)*p - 4*(-k1*k1+k0*k0)*p*p;
    if(top==0 || bot==0){std::cout<<"ERROR LpB: top: "<<top<<" bot: "<<bot<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
	return gsl_sf_log_abs(top/bot);
}

double LmB(double p, double M, double k0, double k1){
	double top = pow(k0*k0-k1*k1-M*M,2)-4*pow(k0+k1,2)*p*p;
	double bot = pow(k0*k0-k1*k1-M*M,2)-4*pow(k0-k1,2)*p*p;
    if(top==0||bot==0){std::cout<<"ERROR LmB: top: "<<top<<" bot: "<<bot<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
	return gsl_sf_log_abs(top/bot);
}

double LpF(double p, double M, double k0, double k1){

        double top=pow(k0*k0-k1*k1+M*M,2)+4*k1*(-k1*k1+k0*k0+M*M)*p-4*((p*p+M*M)*k0*k0-p*p*k1*k1);
        double bot=pow(k0*k0-k1*k1+M*M,2)-4*k1*(-k1*k1+k0*k0+M*M)*p+4*((p*p+M*M)*k0*k0-p*p*k1*k1);
    if(top==0||bot==0){std::cout<<"ERROR LpF: top: "<<top<<" bot: "<<bot<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
         return gsl_sf_log_abs(top/bot);

}

double LmF(double p, double M, double k0, double k1){
    
        double top=pow(k0*k0-k1*k1+M*M,2)-4*(sqrt(p*p+M*M)*k0+p*k1)*(sqrt(p*p+M*M)*k0+p*k1);
        double bot=pow(k0*k0-k1*k1+M*M,2)-4*(sqrt(p*p+M*M)*k0-p*k1)*(sqrt(p*p+M*M)*k0-p*k1);
    if(top==0||bot==0){std::cout<<"ERROR LmF: top: "<<top<<" bot: "<<bot<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
	return gsl_sf_log_abs(top/bot);
}


// why is void here?
double TRKS_integrand(double p, void *d){
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k1 = params->k1;
	double k0 = params->k0;
	double T = params->T;

	//return nB(p,T)/(8*PI*PI)*p*(4+(M*M-k1*k1+k0*k0)*LpB(p,M,k0,k1)/(2*p*k1))+ (nF(p,T,M)*p*p)/(2*2*PI*2*PI*sqrt(p*p+M*M))*(4-(M*M-k1*k1+k0*k0)*LpF(p,M,k0,k1)/(2*p*k1));
    return 1/(8*PI*PI)*(p*nB(p,T)*(4+((k0*k0-k1*k1+M*M)/(2*p*k1))*LpB(p,M,k0,k1))+  ((p*p*nF(p,T,M))/(p*p+M*M))*(4-((k0*k0-k1*k1+M*M)/(2*p*k1))*LmF(p,M,k0,k1)));
}


double TRKS(double T, double M,double k0, double k1)
{

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	struct dispParams params = {T,M,k0,k1};
    double result,error;
	gsl_function F;
	/*
	 *Base is the generic starting value for what I base the integral size off.
	 * Then this increases temperatyre in steps of "T" until I have reached a value below that the predefined threshold, T_TRESH
	 * It then integrates this from 0 to value of T where cross the T_Treshold.
	 */
	double base = TRKS_integrand(T,&params);
	double iT = 2;

	for(iT=2;fabs(TRKS_integrand(iT*T,&params)/base)>=T_TRESHOLD;iT=iT+1)
	{
        //if(DEBUG_MODE){std::cout<<std::setprecision(9)<<" # TRKS integral base "<<base<<" iT: "<<iT<<" rat:  "<<fabs(TRKS_integrand(iT*T,&params)/base)<<std::endl;}
	}


	F.function = &TRKS_integrand;
	F.params = &params;
	
	gsl_integration_qags(&F,0.0, iT*T, ABS, REL, 5000, w, &result, &error);
	gsl_integration_workspace_free (w);
    

	return result;
   

}


double TRUS_integrand(double p, void *d){
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k1 = params->k1;
	double k0 = params->k0;
	double T = params->T;

	//return nB(p,T)/(2*2*PI*2*PI*k1)*(p*LmB(p,M,k0,k1)+k0*LpB(p,M,k0,k1));
    //return 1/(2*2*PI*2*PI)*(  p*LmB(p,M,k0,k1)*nB(p,T)/k1 + k0*LpB(p,M,k0,k1)*nB(p,T)/k1 + p*nF(p,T,M)*LmF(p,M,k0,k1)/k1 );
    return 1/(8*PI*PI)*(     (p*nB(p,T)/k1)* (LmB(p,M,k1,k0)+(k0/p)*LpB(p,M,k1,k0) ) +( p*nF(p,T,M)/k1)*( LmF(p,M,k0,k1)) );
}


double TRUS(double T, double M,double k0, double k1)
{
    //disParams getting defined again, can it be defined globally?
    
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	struct dispParams params = {T,M,k0,k1};
	double result,error;
	gsl_function F;
	F.function = &TRUS_integrand;
	F.params = &params;
	
	double base = TRUS_integrand(T,&params);
	double iT = 2;

	for(iT=2;fabs(TRUS_integrand(iT*T,&params)/base)>=T_TRESHOLD;iT=iT+1)
	  {
	 // if(DEBUG_MODE){	std::cout<<std::setprecision(9)<<"# TRKU integral base "<<base<<" iT: "<<iT<<" rat:  "<<fabs(TRUS_integrand(iT*T,&params)/base)<<std::endl;}
		}
    //std::cout << "now: " << iT << std::endl;
	gsl_integration_qags(&F,0.0, iT*T, ABS, REL, 5000, w, &result, &error);
	gsl_integration_workspace_free (w);

	return result;

}

double TRS_integrand(double p, void*d){
	
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k1 = params->k1;
	double k0 = params->k0;
	double T = params->T;


	return -1.0/(4.0*PI*PI)*p/k1*(1.0/k0*nF(p,T,M)*LpF(p,M,k0,k1)-1.0/p*nB(p,T)*LpB(p,M,k0,k1));


}

double TRS(double T,double M,double k0,double k1)
{

    
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	struct dispParams params = {T,M,k0,k1};
	double result,error;
	gsl_function F;
	F.function = &TRS_integrand;
	F.params = &params;
	
	double base = TRS_integrand(T,&params);
	double iT = 2;

	for(iT=2;fabs(TRS_integrand(iT*T,&params)/base)>=T_TRESHOLD;iT=iT+1)
	  {
	 // if(DEBUG_MODE){	std::cout<<std::setprecision(9)<<"# TRKU integral base "<<base<<" iT: "<<iT<<" rat:  "<<fabs(TRUS_integrand(iT*T,&params)/base)<<std::endl;}
		}
    //std::cout << "now: " << iT << std::endl;
	gsl_integration_qags(&F,0.0, iT*T, ABS, REL, 5000, w, &result, &error);
	gsl_integration_workspace_free (w);

	return result;

}



double sigmaA(double T, double M, double k0, double k1){
	return 1/(k1*k1)*(TRKS(T,M,k0,k1)-k0*TRUS(T,M,k0,k1));
}


double sigmaB(double T, double M, double k0, double k1){
	return (pow(k0/k1,2)-1)*TRUS(T,M,k0,k1)-(k0/pow(k1,2))*TRKS(T,M,k0,k1);
   }

double sigmaC(double T, double M, double k0, double k1){
	return TRS(T,M,k0,k1);
}

//opposite now, put in k0 get k1
double dispEquation(double k0, void * d){
	struct dispParams2 * params = (struct dispParams2 *)d;
	double T = params->T;
	double M = params->M;
        double k1 = params->k1;
    
    
        //if(DEBUG_MODE){std::cout<<"sigma"<<" "<<(1+sigmaA(*params))*k0+sigmaB(*params)+(1+sigmaA(*params))*w1<<std::endl;}
         // double ans_pos = (1+sigmaA(*params))*w0+sigmaB(*params)+(1+sigmaA(*params))*k1;
         // double ans_neg = (1+sigmaA(T,M,k0,k1))*k0+sigmaB(T,M,k0,k1)-(1+sigmaA(T,M,k0,k1))*k1;
	  double ans_mass= (1+sigmaA(T,M,k0,k1))*k0+sigmaB(T,M,k0,k1)-sqrt(pow(k1*(1.0+sigmaA(T,M,k0,k1)),2)+pow(M*(1-sigmaC(T,M,k0,k1)),2) );
   	 return ans_mass;
   // return ans_pos;

}

double dispSolved(double T, double M, double k1){

	int iter = 0, max_iter = 250;
   	int status;
	double r = 0;
	// INitialise gsl root solving algorithm, going to use brent as its a brackeing algorithm that as long as we get a point on either sie, is guarrenteed to
	// find the contained root.
	
	const gsl_root_fsolver_type * S   = gsl_root_fsolver_brent;
	gsl_root_fsolver * s     = gsl_root_fsolver_alloc (S);
    struct dispParams2 p = {T,M,k1};

	//This for loop will start a bit to the left and roight of input k0 and increas:e the bracket until we get opposite signs
	//Will only fail if its all negative or positive, which can happen at very low k/T (presumably due to not complete expressions)
	
	double x_lo=k1;
	double x_hi=k1;
	

	double m = 1e-8;
    for(m = 1e-6; whatsthatsign(dispEquation(k1*pow(10,m),&p))==whatsthatsign(dispEquation(k1*pow(10,-m),&p)); m=m+0.005)
   	{
        if(DEBUG_MODE){std::cout<<m<<" "<<pow(10,m)<<" "<<whatsthatsign(dispEquation(k1*pow(10,m),&p))<<" "<<whatsthatsign(dispEquation(k1*pow(10,-m),&p))<<std::endl;}
        if(DEBUG_MODE){std::cout<<m<<" "<<pow(10,m)<<" "<<dispEquation(k1*pow(10,m),&p)<<" "<<-m<<" "<<pow(10,-m)<<" "<<dispEquation(k1*pow(10,-m),&p)<<std::endl;}
		if(m>4)
		{
			std::cout<<"# ERROR: dispSolved@'dispersion.c' rootfinder: Wandered from k0, didnt bracket the root. "<<std::endl;
		}
	
	}
		//std::cout<<m*multiplier<<" "<<dispEquation(w0*(m*multiplier),&p)<<" "<<dispEquation(w0/(m*multiplier),&p)<<std::endl;
		x_lo=std::min(k1*pow(10,-m),k1*pow(10,m));
		x_hi=std::max(k1*pow(10,-m),k1*pow(10,m));


	//Same as with integration, we need to set it up in terms of a gsl_function F .
	gsl_function F;
	F.function =&dispEquation;
	F.params =&p;

	gsl_root_fsolver_set(s,&F,x_lo,x_hi);
	

	// In general dont really need much iterations, very efficient algorithm
	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_x_lower(s);
		x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo,x_hi,0,1e-8);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
    
   // if(iter>max_iter-2){std::cout<<" Hitting max iter?"<<std::endl;}

	// Always free your fsolver for memory reasons
	gsl_root_fsolver_free (s); 


    if(r<0){std::cout<<"#ERROR: k0 is negative energy state"<<std::endl;}
    
	return r;

}

