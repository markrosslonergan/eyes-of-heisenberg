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
#include "dispersion.h"

int main(int argc, char * argv[])
{

double cutoff = 0.0;
    double M = 1e13;
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


	//std::cout<<"#Cutoff: "<<cutoff<<std::endl;
	double baseT = 1e15;

if(boolMe){
     //T,M,k1,k0
	

	//for(double x = -2; x<1; x=x+0.0033)
    for(double x = 0.6;x >0.04 ; x=x-0.01)
	{
        double k1 = x*baseT;
        // dispSolved now takes (T,M,k1) simple.
        double k0_09 = dispSolved(baseT, 0.9*baseT, k1);
        std::cout<<fabs(k1)/baseT<<" "<<(k0_09-fabs(k1))/baseT<<" "<<x<<std::endl;
	}


} else {

	//for(double z = -2; z<1; z=z+0.05){
	//	std::cout<<pow(10,z)<<" "<<epsilon_eff_ratio(M/(pow(10,z)),M,cutoff)-1.0<<std::endl;
	//}
    double ik0 = 1.0;
    baseT=1e14;
    struct dispParams2 mypar ={baseT,0.9*baseT,0.4*baseT};
    for(ik0=-2.0; ik0<0.0;ik0=ik0+0.01){
        std::cout<<pow(10,ik0)<<"  "<<dispEquation(pow(10,ik0)*baseT,&mypar)/baseT<<std::endl;
        
        }
    }
    
    

return 0;
}
