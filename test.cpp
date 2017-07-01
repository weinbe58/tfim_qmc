

#include "qaqmc.h"
#include <iostream>



int main(){
	// cast simulation object:
	const int L=5;
	const int d=1;
	const int M=32;
	const int Nm=0;
	const int mstep=1000000;
	const double S_i=1.0;
	const double S_f=0.0;
	std::string bc = "pzpz";
	std::string file = "output.dat";

	qaqmc<L,d,M,Nm,mstep> qmc(S_i,S_f,bc,file);

	for(int i=0;i<3;i++){
		std::cout << "eq:  " << i+1 << std::endl;
		qmc.EQstep();
	}

	const int nbin = 1000000;
	for(int i=0;i<nbin;i++){
		std::cout << "bin:  " << i+1 << std::endl;
		qmc.MCstep();
		// qmc.print_opstr(true);
		qmc.write_out();
	}

	return 0;

}








