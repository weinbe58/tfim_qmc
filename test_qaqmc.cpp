

#include "qaqmc.h"
#include <iostream>



int main(){
	// cast simulation object:
	const int L=10;
	const int d=1;
	const int M=10000;
	const int Nm=100;
	const int mstep=10000;
	const double S_i=1.0;
	const double S_f=0.0;
	std::string bc = "pzpz";
	std::string file = "test_output.dat";

	qaqmc<L,d,M,Nm,mstep> q(S_i,S_f,bc,file);

	for(int i=0;i<3;i++){
		std::cout << "eq:  " << i+1 << std::endl;
		q.EQstep();
	}

	const int nbin = 100;
	for(int i=0;i<nbin;i++){
		std::cout << "bin:  " << i+1 << std::endl;
		q.MCstep();
		// q.print_opstr(true);
		q.write_out();
	}

	return 0;

}








