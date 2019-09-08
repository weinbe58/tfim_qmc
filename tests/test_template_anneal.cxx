

#include "template_anneal.h"
#include <iostream>



int main(){
	// cast simulation object:
	const int L=100;
	const int d=1;
	const int N = get_N<L,d>::N;
	const int M=1000;
	const int mann=1000;
	const double S_f=0.5;
	std::string file = "test_output.dat";

	qaqmc<L,d,M,mann> q(S_f,file);

	const int nbin = 100;
	for(int i=0;i<nbin;i++){
		std::cout << "bin:  " << i+1 << std::endl;
		q.MCramp(1000);
		// q.print_opstr(true);
		q.write_out();
	}

	return 0;

}








