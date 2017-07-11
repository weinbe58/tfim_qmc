

#include "qmc.h"
#include <iostream>



int main(){
	// cast simulation object:
	const int L=100;
	const int d=1;
	const int M=1000;
	const int Nm=100;
	const int mstep=10000;
	const double S=0.4;
	std::string bc = "pzpz";
	std::string file = "test_output.dat";

	qmc<L,d,M,Nm,mstep> q(S,bc,file);

	for(int i=0;i<3;i++){
		std::cout << "eq:  " << i+1 << std::endl;
		q.EQstep();
	}

	const int nbin = 100;
	for(int i=0;i<nbin;i++){
		std::cout << "bin:  " << i+1 << std::endl;
		q.MCstep();
		// q.print_opstr(true);
		// q.write_out();
	}

	return 0;

}








