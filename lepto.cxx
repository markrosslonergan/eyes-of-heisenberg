#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "epsilon.h"
#include "washout.h"

int main(int argc, char * argv[])
{

double cutoff = 0.0;
double M = 1e11;
bool boolMe = false;

	int c;
	while ((c = getopt (argc, argv, "WC:M:")) != -1)
   	{ 
	switch (c) 
      	{
      		case 'C':
			cutoff =strtof(optarg, NULL); 
        		break;
      		case 'M':
			M=strtof(optarg,NULL);
			break;
		case 'W':
			boolMe=true;
			break;
      		
      		case '?':
			printf("Allowed arguments:\n"
				"\t-C\tGive the lower bound cutoff for integral in epsilon_eff, Default 0.0.\n"		
				"\t-M\tInput sterile Mass in GeV, default 1e11.\n"
				"\t-W\tWashout run, default false, toggles on.\n"
				);
          	return 1;
      		default:
			printf("I don't know how you got here.\n");
        		abort ();
 	}}


	std::cout<<"#Cutoff: "<<cutoff<<std::endl;


if(boolMe){
	
	/*	double T = 1e12;
	struct washoutParams params = {T,M};
	for(double x = 0; x<20; x=x+0.2)
	{
		std::cout<<x<<" "<<washout_integrand(x,&params)<<std::endl;
	}*/
	for(double z = -2; z<1; z=z+0.05){
		std::cout<<pow(10,z)<<" "<<washout(1,M/(pow(10,z)),M)<<std::endl;
	}

} else {
	/*double T = 1e12;
	struct epsilonParams params = {T,M};
	for(double kk = 4; kk<13.5; kk=kk+0.1)
	{
		std::cout<<pow(10,kk)<<" "<<epsilon_numerator_integrand (pow(10,kk),&params)<<" "<<epsilon_denominator_integrand(pow(10,kk),&params)<<std::endl;
	}*/

	for(double z = -2; z<1; z=z+0.05){
		std::cout<<pow(10,z)<<" "<<epsilon_eff_ratio(M/(pow(10,z)),M,cutoff)-1.0<<std::endl;
	}
}
return 0;
}
