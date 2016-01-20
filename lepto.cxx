#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "epsilon.h"
#include "washout.h"
#include "dispersion.h"

int main(int argc, char * argv[])
{


/******************************************************************
 *
 *	Inputed Parameters and Boolian quantaties, with defaults!
 *
 ******************************************************************
 */

double M = 2;
double T = 1;
bool boolDisp = false;
bool use_nf = false;
bool use_c = false;
bool test_mode = false;
/******************************************************************
 *
 *	Main Flag seperation switch statement, run -h for options
 *
 ******************************************************************
 */

	int c;
	while ((c = getopt (argc, argv, "DM:T:-:")) != -1)
   	{ 
	switch (c) 
      	{
		case 'M':
			M = strtof(optarg,NULL);
			break;
		case 'T':
			T = strtof(optarg,NULL);
			break;
		case 'D':
			boolDisp=true;
			break;
		case '-':
			if(!strcmp(optarg,"dispersion"))
			{
				boolDisp = true;
			}else if(!strcmp(optarg,"using-nf"))
			{
				use_nf = true;
			}else if(!strcmp(optarg,"using-c"))
			{
				use_c = true;
			}else if(!strcmp(optarg,"test"))
			{
				test_mode = true;
			}else
			{
				std::cout<<"Passed a -- flag that shouldnt be!"<<std::endl;
				return -1;
			}
			break;
      		
      		case '?':
			printf("Allowed arguments:\n"
				"\t-M\tInput sterile Mass in GeV, default 2.\n"
				"\t-T\tInput temperature in GeV, default 1.\n"
				"\t-D --dispersion\t Run a dispersion plot, for inputted M and T.\n"
				"\t--using-nf\t Use fermionic Nf terms in LpmB, default off.\n"
				"\t--using-c\t Use mass dependant sigma C terms in dispersion eq, default off\n"
				"\t--test\t run a standard MoT 2,1.5,1.2,1,0.9,0.7 dispersion.\n"
				);
          	return 1;
      		default:
			printf("I don't know how you got here.\n");
        		abort ();
 	}}

/******************************************************************
 *
 *	Derived Parameters and Boolian quantaties
 *
 ******************************************************************
 */
	double 	mot = M/T;


/******************************************************************
 *
 *	Main function decision tree
 *
 ******************************************************************
 */
if(test_mode){

	std::vector<double > list_mot = {2,1.5,1.2,1,0.9,0.7};
	
	for(int i =0; i< list_mot.size();i++){

		double iMoT=list_mot[i];
		double iT=1;
		double iM = iMoT*iT;

		std::stringstream str; str.str(std::string()); str.clear();
       		std::ofstream myfile;
		
		//Filename in tests/dispersion_test
		str<<std::fixed<<std::setprecision(2)<<"tests/dispersion_test/MoT_"<<iM/iT;
		if(use_c){
			str<<"_wC";
		}
		if(use_nf){
			str<<"_wNF";
			
		}
		str<<".dat";

		myfile.open(str.str().c_str());
	
		std::cout<<"Starting M/T of : "<<iM/iT<<" to tests/dispersion_test/*.dat"<<std::endl;
	
        	for(double k1oT = 1.0;k1oT > 0.02 ; k1oT=k1oT-0.0101)
		{	
			double k1 = k1oT*iT;
		        double found_k0 = dispSolved(iT, iM, k1, use_c, use_nf);

        		myfile<<fabs(k1oT)<<" "<<(found_k0-fabs(k1))/iT<<std::endl;
		}


		myfile.close();
	}

		std::cout<<"Running gnuplot script to plot tests/dispersion_test/dispersion_plot.png"<<std::endl;
		system("./tests/dispersion_test/plot_dispersion.sh");
		std::cout<<"Done!"<<std::endl;


}else if(boolDisp){


	/*if(false){ // plotting integrand

		//T M k1 k0
		double Mm = 1.2*baseT;
		double k11 = 0.3*baseT;
		double k00 = 0.04*baseT+k11;
		struct dispParams my = {baseT, Mm,k00,k11,use_nf};

		for(double p = 12; p<=16; p=p+0.00103120){

			std::cout<<pow(10,p)<<" "<<TRKS_integrand(pow(10,p),&my)<<" "<<nB(pow(10,p),baseT)<<"  "<<TRUS_integrand(pow(10,p),&my)<<" "<<LpB_negative(pow(10,p),Mm,k00,k11)<<std::endl;
		}
	}*/
	double kmin=0.02;
	double kmax=1;
	plot_disp(T,M,use_c,use_nf,kmin,kmax);


} else {



	//for(double z = -2; z<1; z=z+0.05){
	//	std::cout<<pow(10,z)<<" "<<epsilon_eff_ratio(M/(pow(10,z)),M,cutoff)-1.0<<std::endl;
	//}
    double ik0 = 1.0;
double    baseT=1.0;
    double MoT = 0.9;
    struct dispParams2 mypar ={baseT,MoT*baseT,0.05*baseT,use_c, use_nf};
    struct dispParams2 mypar2 ={baseT,MoT*baseT,0.1*baseT,use_c, use_nf};
    struct dispParams2 mypar2b={baseT,MoT*baseT,0.13*baseT,use_c, use_nf};
    struct dispParams2 mypar3 ={baseT,MoT*baseT,0.15*baseT,use_c,use_nf};
    struct dispParams2 mypar4 ={baseT,MoT*baseT,0.2*baseT,use_c, use_nf};


    for(ik0=0.8; ik0>0.01;ik0=ik0-0.001002103){

    struct dispParams myparNEW = {baseT,MoT*baseT,ik0,0.1, use_nf};
	//Dispersion Curves
    //    std::cout<<pow(10,ik0)<<"  "<<dispEquation(pow(10,ik0)*baseT,&mypar)/baseT<<" "<<dispEquation(pow(10,ik0)*baseT,&mypar2)/baseT<<" "<<dispEquation(pow(10,ik0)*baseT,&mypar2b)/baseT<<" "<<dispEquation(pow(10,ik0)*baseT,&mypar3)/baseT<<" "<<dispEquation(pow(10,ik0)*baseT,&mypar4)/baseT<<std::endl;

//    std::cout<<ik0<<" "<<dispEquation(ik0,&mypar2)<<" "<<TRKS_integrand(0.5,&myparNEW)<<"  "<<TRUS_integrand(0.5,&myparNEW)<<" "<<dispEquationN(ik0,&mypar2)<<" "<<sigmaA(1,0.9,ik0,0.1)<<" "<<sigmaB(1,0.9, ik0,0.1)<<std::endl;    
  
  
   //// Sigma A
       //
        // std::cout<<pow(10,ik0)<<"  "<<sigmaA(baseT,mypar.M,pow(10,ik0)*baseT,mypar.k1)<<" "<<sigmaA(baseT,mypar2.M,pow(10,ik0)*baseT,mypar2.k1)<<" "<<sigmaA(baseT,mypar2b.M,pow(10,ik0)*baseT,mypar2b.k1)<<" "<<sigmaA(baseT,mypar3.M,pow(10,ik0)*baseT,mypar3.k1)<<" "<<sigmaA(baseT,mypar4.M,pow(10,ik0)*baseT,mypar4.k1)<<std::endl;
	//  sigma B
      // std::cout<<pow(10,ik0)<<"  "<<sigmaB(baseT,mypar.M,pow(10,ik0)*baseT,mypar.k1)<<" "<<sigmaB(baseT,mypar2.M,pow(10,ik0)*baseT,mypar2.k1)<<" "<<sigmaB(baseT,mypar2b.M,pow(10,ik0)*baseT,mypar2b.k1)<<" "<<sigmaB(baseT,mypar3.M,pow(10,ik0)*baseT,mypar3.k1)<<" "<<sigmaB(baseT,mypar4.M,pow(10,ik0)*baseT,mypar4.k1)<<std::endl;
	//TRUS	
      // std::cout<<pow(10,ik0)<<"  "<<TRUS(baseT,mypar.M,pow(10,ik0)*baseT,mypar.k1)<<" "<<TRUS(baseT,mypar2.M,pow(10,ik0)*baseT,mypar2.k1)<<" "<<TRUS(baseT,mypar2b.M,pow(10,ik0)*baseT,mypar2b.k1)<<" "<<TRUS(baseT,mypar3.M,pow(10,ik0)*baseT,mypar3.k1)<<" "<<TRUS(baseT,mypar4.M,pow(10,ik0)*baseT,mypar4.k1)<<std::endl;
	
	//TRKS	
      // std::cout<<pow(10,ik0)<<"  "<<TRKS(baseT,mypar.M,pow(10,ik0)*baseT,mypar.k1)<<" "<<TRKS(baseT,mypar2.M,pow(10,ik0)*baseT,mypar2.k1)<<" "<<TRKS(baseT,mypar2b.M,pow(10,ik0)*baseT,mypar2b.k1)<<" "<<TRKS(baseT,mypar3.M,pow(10,ik0)*baseT,mypar3.k1)<<" "<<TRKS(baseT,mypar4.M,pow(10,ik0)*baseT,mypar4.k1)<<std::endl;

        }
    }
    
    

return 0;
}
