

#include "template_proj.h"
#include <iostream>



int main(){
	// cast simulation object:
	const int L=4;
	const int d=1;
	const int Nm=1;
	const int M=40*get_N<L,d>::N*get_N<L,d>::N+Nm/2;
	const int mstep=10000;
	const double S=0.5;
	std::string bc = "pzpz";
	std::string file = "test_output_2.dat";

	proj<L,d,M,Nm,mstep> q(S,bc,file);

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








