

#include "sse_local.h"
#include <iostream>
#include <iomanip>
#include <tuple>



void sq_lattice(int L,int d,std::vector<int> &bst){
	int N = int(std::pow(L,d));
	for(int i=0;i<N;i++){
		int LL = 1;
		for(int j=0;j<d;j++){
			int id=(i/LL)%L;
			int b = d*i+j;

			bst.push_back(i);
			bst.push_back(i+((id+1)%L-id)*LL);
			LL *= L;
		}
	}
}

void EQsweep(sse_local &qmc,int mstep){
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		qmc.check_M();
	}
}

void MCsweep(sse_local &qmc,int mstep,double &ma,double &m2,double &m4){
	ma = m2 = m4 = 0;
	int N = qmc.get_N();
	int Nb = qmc.get_Nb();
	int M = qmc.get_M();
	std::vector<signed char> spins(N,0);
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		// qmc.print_opstr(false);
		for(auto tup=qmc.state_begin();tup!=qmc.state_end();tup++){
			double mm = std::get<2>(*tup);
			int n = std::get<1>(*tup);
			ma += n*std::abs(mm);
			m2 += n*std::pow(mm,2);
			m4 += n*std::pow(mm,4);
		}
	}

	ma /= (M*size_t(mstep)+1.1e-15);
	m2 /= (M*size_t(mstep)+1.1e-15);
	m4 /= (M*size_t(mstep)+1.1e-15);
}




int main(int argc, char const *argv[])
{
	int L = 4;
	int d = 1;
	double beta = 0.1;
	double S = 0.0;
	int N = int(std::pow(double(L),d));

	std::vector<signed char> bc;

	for(int i=0;i<N;i++){
		bc.push_back(1);
	}

	std::vector<int> bst;

	sq_lattice(L,d,bst);

	double eJ,ma,m2,m4;

	// proj_local qmc(2*M,N,d*N,bst.data(),S,1,1,bc,bc);
	sse_local qmc(beta,N,d*N,bst.data(),S);
	for(int i=0;i<5;i++){
		EQsweep(qmc,10000);
	}

	for(int i=0;i<10;i++){
		MCsweep(qmc,100000,ma,m2,m4);
		std::cout << std::scientific << std::setprecision(10) << std::setw(20);
		std::cout << ma/N << std::setw(20);
		std::cout << m2/(N*N) << std::setw(20);
		std::cout << m4/(N*N*N*N) << std::setw(20);
		std::cout << std::endl;
	}


	return 0;
}


