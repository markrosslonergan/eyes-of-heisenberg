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
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/multiprecision/cpp_dec_float.hpp>


using namespace boost::multiprecision;

#define GW 2.0
#define PI 3.14149
//#define REL 0.000001//0.0001
#define REL 1e-7
#define ABS 0.0
//#define ABS 0
#define DEBUG_MODE false
//#define T_TRESHOLD 1e-10
#define T_TRESHOLD 1e-12

#define P_MIN 0.0
#define T_MAX iT*T 

//#define USE_C false Depreciated code, now flagged
//#define USE_NF 0.0 Depreciated code
#define NUM_INT_SPACE 50000


double nB(double p, double T){
    double ans = 0.0;
//   if (p < T_MAX){ ans = 1.0/(exp(p/T)-1.0);}
	 ans = 1.0/(exp(p/T)-1.0); 
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


/*
double LpB(double p, double M, double k0, double k1){

	double K = 

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
	        double bot=pow(k0*k0-k1*k1+M*M,2)-4*k1*(-k1*k1+k0*k0+M*M)*p-4*((p*p+M*M)*k0*k0-p*p*k1*k1);
	    if(top==0||bot==0){std::cout<<"ERROR LpF: top: "<<top<<" bot: "<<bot<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
           return gsl_sf_log_abs(top/bot);

}

double LmF(double p, double M, double k0, double k1){
	    
	        double top=pow(k0*k0-k1*k1+M*M,2)-4*(sqrt(p*p+M*M)*k0+p*k1)*(sqrt(p*p+M*M)*k0+p*k1);
		double bot=pow(k0*k0-k1*k1+M*M,2)-4*(sqrt(p*p+M*M)*k0-p*k1)*(sqrt(p*p+M*M)*k0-p*k1);
		if(top==0||bot==0){std::cout<<"ERROR LmF: top: "<<top<<" bot: "<<bot<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
	return gsl_sf_log_abs(top/bot);
}
*/

double LpB_pole_top(double M, double k0, double k1){

	return fabs((k1*k1-k0*k0-M*M)/(2*k0+2*k1));
}
double LpB_pole_bot(double M, double k0, double k1){
	//if(k0==k1){std::cout<<"# ERROR: LpB(M,ko,k1) division by zero"<<std::endl;}
	return fabs((k1*k1-k0*k0-M*M)/(2*k0-2*k1));
}




double LpB(double p, double M, double k0, double k1){
	double KK = k0*k0-k1*k1;
	
	double top1 = KK-M*M+2*p*(k0+k1);
	double bot1 = KK-M*M+2*p*(k0-k1);

	double top2 = KK-M*M-2*p*(k0-k1);
	double bot2 = KK-M*M-2*p*(k0+k1);
  	
	double poleT = fabs(LpB_pole_top(M,k0,k1));
	double poleB = fabs(LpB_pole_bot(M,k0,k1));

	if (p==poleT || p==poleB){std::cout<<"# ERROR: I am directly calling a pole in TRKS_integrand"<<std::endl;}
	
	//if(top1*top2 == 0 || bot1*bot2==0){std::cout<<"ERROR LpB: top: "<<top1*top2<<" bot: "<<bot1*bot2<<" k0: "<<k0<<" k1: "<<k1<<std::endl;
	//};

//	if(top1*top2/(bot1*bot2)<0){ std::cout<<"ERROR LpB: negative log.  top1: "<<top1<<" top2: "<<top2<<" bot1: "<<bot1<<" bot2 : "<<bot2<<" k0: "<<k0<<" k1: "<<k1<<" p: "<<p<<" M: "<<M<<std::endl;}

//	};
	//if((whatsthatsign(top1/bot1)!=whatsthatsign(top2/bot2)) ){ std::cout<<"ERROR LpB: mismatched complex.  top1: "<<top1<<" top2: "<<top2<<" bot1: "<<bot1<<" bot2 : "<<bot2<<" k0: "<<k0<<" k1: "<<k1<<std::endl;
	//};


//	return gsl_sf_log_abs(;
	return gsl_sf_log_abs((top1*top2)/(bot1*bot2));
}

double LpB_negative(double p, double M, double k0, double k1){
	double KK = k0*k0-k1*k1;
	
	double top1 = KK-M*M+2*p*(k0+k1);
	double bot1 = KK-M*M+2*p*(k0-k1);

	double top2 = KK-M*M-2*p*(k0-k1);
	double bot2 = KK-M*M-2*p*(k0+k1);
  	
	double poleT = fabs(LpB_pole_top(M,k0,k1));
	double poleB = fabs(LpB_pole_bot(M,k0,k1));

//	if(top1*top2/(bot1*bot2)<0){ std::cout<<"ERROR LpB2: negative log.  top1: "<<top1<<" top2: "<<top2<<" bot1: "<<bot1<<" bot2 : "<<bot2<<" k0: "<<k0<<" k1: "<<k1<<" p: "<<p<<" M: "<<M<<std::endl;}

//	if (p==poleT || p==poleB){std::cout<<"# ERROR: I am directly calling a pole in TRKS_integrand"<<std::endl;}
	
	//if(top1*top2 == 0 || bot1*bot2==0){std::cout<<"ERROR LpB: top: "<<top1*top2<<" bot: "<<bot1*bot2<<" k0: "<<k0<<" k1: "<<k1<<std::endl;
	//};

//	if(top1*top2/(bot1*bot2)<0){ std::cout<<"ERROR LpB: negative log.  top1: "<<top1<<" top2: "<<top2<<" bot1: "<<bot1<<" bot2 : "<<bot2<<" k0: "<<k0<<" k1: "<<k1<<" p: "<<p<<" M: "<<M<<std::endl;}
//	};
	//if((whatsthatsign(top1/bot1)!=whatsthatsign(top2/bot2)) ){ std::cout<<"ERROR LpB: mismatched complex.  top1: "<<top1<<" top2: "<<top2<<" bot1: "<<bot1<<" bot2 : "<<bot2<<" k0: "<<k0<<" k1: "<<k1<<std::endl;
	//};


//	return gsl_sf_log_abs(;
	return whatsthatsign((top1*top2)/(bot1*bot2));
	//return 323.0;
}

double LmB(double p, double M, double k0, double k1){

	double KK = k0*k0-k1*k1;
	
	double top1 = KK-M*M+2*p*(k0+k1);
	double bot1 = KK-M*M+2*p*(k0-k1);

	double top2 = KK-M*M-2*p*(k0-k1);
	double bot2 = KK-M*M-2*p*(k0+k1);

	double poleT = fabs(LpB_pole_top(M,k0,k1));
	double poleB = fabs(LpB_pole_bot(M,k0,k1));

	if (p==poleT || p==poleB){std::cout<<"# ERROR: I am directly calling a pole in TRKS_integrand"<<std::endl;}
	
	//if((whatsthatsign(top1)!=whatsthatsign(bot1)) && (whatsthatsign(top2)!=whatsthatsign(bot2) ) ){ std::cout<<"ERROR LmB: mismatched complex.  top1: "<<top1<<" top2: "<<top2<<" bot1: "<<bot1<<" bot2 : "<<bot2<<" k0: "<<k0<<" k1: "<<k1<<std::endl;
	//};


//	if(top1*bot2 == 0 || bot1*top2==0){std::cout<<"ERROR LmB: top: "<<top1*bot2<<" bot: "<<bot1*top2<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
	return gsl_sf_log_abs((top1*bot2)/(bot1*top2));

}

double LpF(double p, double M, double k0, double k1){

	double KK = k0*k0-k1*k1;
	double E = sqrt(p*p+M*M);
	
	double top1 = KK+M*M+2*E*k0+2*p*k1;
	double bot1 = KK+M*M+2*E*k0-2*p*k1;

	double top2 = KK+M*M-2*E*k0+2*p*k1;
	double bot2 = KK+M*M-2*E*k0-2*p*k1;
	
	return gsl_sf_log_abs(top1*top2/(bot1*bot2));

    //if(top==0||bot==0){std::cout<<"ERROR LpF: top: "<<top<<" bot: "<<bot<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};

}

double LmF(double p, double M, double k0, double k1){
 
	double KK = k0*k0-k1*k1;
	double E = sqrt(p*p+M*M);
	
	double top1 = KK+M*M+2*E*k0+2*p*k1;
	double bot1 = KK+M*M+2*E*k0-2*p*k1;

	double top2 = KK+M*M-2*E*k0+2*p*k1;
	double bot2 = KK+M*M-2*E*k0-2*p*k1;
	
//	std::cout<<" top "<<top1*bot2<<" bot "<<bot1*top2<<" arg "<<top1*bot2/(bot1*top2)<<std::endl;	
	return gsl_sf_log_abs(top1*bot2/(bot1*top2));

   
    //if(top==0||bot==0){std::cout<<"ERROR LmF: top: "<<top<<" bot: "<<bot<<" k0: "<<k0<<" k1: "<<k1<<std::endl;};
}

double my_log_test(){
	
	double k0 = 2.23872e14;
	double k1 = 9.88551e8;
	double p = 4.91323e15;
	double M = 1.5e15;
	double T = 1e15;

	double mytest = LmF(p,M,k0,k1);
	std::cout<<mytest<<std::endl;

	 //Operations at fixed precision and full numeric_limits support:
	//	 cpp_dec_float_100 b = 2;
	 //std::cout << std::numeric_limits<cpp_dec_float_100>::digits << std::endl;
	 // Note that digits10 is the same as digits, since we're base 10! :
	 // std::cout << std::numeric_limits<cpp_dec_float_100>::digits10 << std::endl;
	  // We can use any C++ std lib function, lets print all the digits as well:
	 // std::cout << std::setprecision(std::numeric_limits<cpp_dec_float_100>::max_digits10)<< log(b) << std::endl; // print log(2)
	 //    // We can also use any function from Boost.Math:
	 //    // These even work when the argument is an expression template:
return 0;
}

double TRKS_integrand(double p, void *d){
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k1 = params->k1;
	double k0 = params->k0;
	double T = params->T;

	bool use_nf = params->use_nf;

	double KK = k0*k0-k1*k1;
	double E = sqrt(p*p+M*M);



	double pre = 0.5*p*p/(2*2*PI*PI);

	double Bbit = (4.0+(KK+M*M)/(2*p*k1)*LpB(p,M,k0,k1))*nB(p,T)/p;
	double Fbit =0;
        if(use_nf){
		Fbit=	(4.0-(KK+M*M)/(2*p*k1)*LpF(p,M,k0,k1))*nF(p,T,M)/E;
	}

	return pre*(Bbit+Fbit);

}


double TRKS(double T, double M,double k0, double k1,bool use_nf)
{

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (NUM_INT_SPACE );
	struct dispParams params = {T,M,k0,k1,use_nf};
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
	
	gsl_integration_qags(&F,P_MIN, T_MAX, ABS, REL, NUM_INT_SPACE, w, &result, &error);
	gsl_integration_workspace_free (w);
    

	return result;
   

}


double TRUS_integrand(double p, void *d){
	struct dispParams * params = (struct dispParams *)d;
	double T = params->T;
	double M = params->M;
	double k0 = params->k0;
	double k1 = params->k1;

	bool use_nf = params->use_nf;
	double KK = k0*k0-k1*k1;
	double E = sqrt(p*p+M*M);

	double preB = 0.5*p/(2*2*PI*PI*k1);

	double Bbit=(LmB(p,M,k0,k1)+k0/p*LpB(p,M,k0,k1))*nB(p,T);
	double Fbit =0;
        if(use_nf){
	 	Fbit=LmF(p,M,k0,k1)*nF(p,T,M);
	}

    	return preB*(Bbit+Fbit);

	//return nB(p,T)/(2*2*PI*2*PI*k1)*(p*LmB(p,M,k0,k1)+k0*LpB(p,M,k0,k1));
    //return 1/(2*2*PI*2*PI)*(  p*LmB(p,M,k0,k1)*nB(p,T)/k1 + k0*LpB(p,M,k0,k1)*nB(p,T)/k1 + p*nF(p,T,M)*LmF(p,M,k0,k1)/k1 );
}





double TRUS(double T, double M,double k0, double k1,bool use_nf)
{
    //disParams getting defined again, can it be defined globally?
    
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (NUM_INT_SPACE);
	struct dispParams params = {T,M,k0,k1,use_nf};
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
	gsl_integration_qags(&F,P_MIN, T_MAX, ABS, REL, NUM_INT_SPACE, w, &result, &error);
	gsl_integration_workspace_free (w);

	return result;

}

double TRS_integrand(double p, void*d){
	
	struct dispParams * params = (struct dispParams *)d;
	double M = params->M;
	double k0 = params->k0;
	double k1 = params->k1;
	double T = params->T;

	bool use_nf = params->use_nf;
	
	double E= sqrt(p*p+M*M);

	double pre = M*p/(2*2*PI*PI*k1);
	double Bbit = -1/p*nB(p,T)*LpB(p,M,k0,k1);
	double Fbit =0;
        
	if(use_nf){
	 	Fbit = 1/E*nF(p,T,M)*LpF(p,M,k0,k1);
	}
	return pre*(Fbit+Bbit);

}

double TRS(double T,double M,double k0,double k1, bool use_nf)
{

    
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (NUM_INT_SPACE);
	struct dispParams params = {T,M,k0,k1,use_nf};
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
	gsl_integration_qags(&F,P_MIN, T_MAX, ABS, REL, NUM_INT_SPACE, w, &result, &error);
	gsl_integration_workspace_free (w);

	return result;

}



double sigmaA(double T, double M, double k0, double k1,bool use_nf){
	return 1/(k1*k1)*(TRKS(T,M,k0,k1,use_nf)-k0*TRUS(T,M,k0,k1,use_nf));
}

double sigmaB(double T, double M, double k0, double k1, bool use_nf){
	return (pow(k0/k1,2)-1)*TRUS(T,M,k0,k1,use_nf)-(k0/pow(k1,2))*TRKS(T,M,k0,k1,use_nf);
   }

double sigmaC(double T, double M, double k0, double k1, bool use_nf){
	return -1/(4*M)*TRS(T,M,k0,k1,use_nf);
}

//opposite now, put in k0 get k1
double dispEquation(double k0, void * d){
	struct dispParams2 * params = (struct dispParams2 *)d;
	double T = params->T;
	double M = params->M;
        double k1 = params->k1;
      	bool use_c = params->use_c;
	bool use_nf= params->use_nf;

	double sC = sigmaC(T,M,k0,k1, use_nf);

	if(!use_c){
		sC=1;
	}

   	double sA = sigmaA(T,M,k0,k1,use_nf);
	double sB = sigmaB(T,M,k0,k1,use_nf);
        //if(DEBUG_MODE){std::cout<<"sigma"<<" "<<(1+sigmaA(*params))*k0+sigmaB(*params)+(1+sigmaA(*params))*w1<<std::endl;}
         // double ans_pos = (1+sigmaA(*params))*w0+sigmaB(*params)+(1+sigmaA(*params))*k1;
         // double ans_neg = (1+sigmaA(T,M,k0,k1))*k0+sigmaB(T,M,k0,k1)-(1+sigmaA(T,M,k0,k1))*k1;
	 // double ans_mass = (1+sigmaA(T,M,k0,k1))*k0+sigmaB(T,M,k0,k1)-sqrt(pow(k1*(1.0+sigmaA(T,M,k0,k1)),2)+pow(M*(1-sigmaC(T,M,k0,k1)),2) );
	 double ans = (1+sA)*k0+sB-sqrt(pow(k1*(1.0+sA),2)+pow(M*(1-sC),2));

   	 return ans;
   // return ans_pos;

}

double dispEquationN(double k0, void * d){
	struct dispParams2 * params = (struct dispParams2 *)d;
	double T = params->T;
	double M = params->M;
        double k1 = params->k1;
     	bool use_c = params->use_c;
	bool use_nf= params->use_nf;

	double sC = sigmaC(T,M,k0,k1, use_nf);

	if(!use_c){
		sC=1;
	}

   	double sA = sigmaA(T,M,k0,k1,use_nf);
	double sB = sigmaB(T,M,k0,k1,use_nf);
    
        //if(DEBUG_MODE){std::cout<<"sigma"<<" "<<(1+sigmaA(*params))*k0+sigmaB(*params)+(1+sigmaA(*params))*w1<<std::endl;}
         // double ans_pos = (1+sigmaA(*params))*w0+sigmaB(*params)+(1+sigmaA(*params))*k1;
         // double ans_neg = (1+sigmaA(T,M,k0,k1))*k0+sigmaB(T,M,k0,k1)-(1+sigmaA(T,M,k0,k1))*k1;
	 // double ans_mass = (1+sigmaA(T,M,k0,k1))*k0+sigmaB(T,M,k0,k1)-sqrt(pow(k1*(1.0+sigmaA(T,M,k0,k1)),2)+pow(M*(1-sigmaC(T,M,k0,k1)),2) );
	double ans = (1+sA)*k0+sB+sqrt(pow(k1*(1.0+sA),2)+pow(M*(1-sC),2));


   	 return ans;
   // return ans_pos;

}

double dispSolved(double T, double M, double k1, double init, bool use_c, bool use_nf){

	int iter = 0, max_iter = 250;
   	int status;
	double r = 0;
	// INitialise gsl root solving algorithm, going to use brent as its a brackeing algorithm that as long as we get a point on either sie, is guarrenteed to
	// find the contained root.
	
	const gsl_root_fsolver_type * S   = gsl_root_fsolver_brent;
	gsl_root_fsolver * s     = gsl_root_fsolver_alloc (S);
        struct dispParams2 p = {T,M,k1,use_c,use_nf};

	//This for loop will start a bit to the left and roight of input k0 and increas:e the bracket until we get opposite signs
	//Will only fail if its all negative or positive, which can happen at very low k/T (presumably due to not complete expressions)
	
	double x_lo=init;
	double x_hi=init;
	

	double m = 1e-8;
    for(m = 1e-6; whatsthatsign(dispEquation(init*pow(10,m),&p))==whatsthatsign(dispEquation(init*pow(10,-m),&p)); m=m+0.0005)
   	{
  //        if(DEBUG_MODE){std::cout<<m<<" "<<pow(10,m)<<" "<<whatsthatsign(dispEquation(k1*pow(10,m),&p))<<" "<<whatsthatsign(dispEquation(k1*pow(10,-m),&p))<<std::endl;}
  //      if(DEBUG_MODE){std::cout<<m<<" "<<pow(10,m)<<" "<<dispEquation(k1*pow(10,m),&p)<<" "<<-m<<" "<<pow(10,-m)<<" "<<dispEquation(k1*pow(10,-m),&p)<<std::endl;}
		if(m>4)
		{
		//	std::cout<<"# ERROR: dispSolved@'dispersion.c' rootfinder: Wandered from k0, didnt bracket the root. "<<std::endl;
		}
	
	}
		//std::cout<<m*multiplier<<" "<<dispEquation(w0*(m*multiplier),&p)<<" "<<dispEquation(w0/(m*multiplier),&p)<<std::endl;
		x_lo=std::min(init*pow(10,-m),init*pow(10,m));
		x_hi=std::max(init*pow(10,-m),init*pow(10,m));


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


double plot_disp(double T, double M, bool use_c, bool use_nf, double k1oT_min, double k1oT_max){
	double init = k1oT_max*T;

        for(double k1oT = k1oT_max;k1oT > k1oT_min ; k1oT=k1oT-0.0101)
	{	
		
		double k1 = k1oT*T;
	        double found_k0 = dispSolved(T, M, k1,init, use_c, use_nf);

        	std::cout<<fabs(k1oT)<<" "<<(found_k0-fabs(k1))/T<<std::endl;
		init = found_k0;
	}


	return 1;
}
