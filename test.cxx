

#include "proj_local.h"
#include "qaqmc_local.h"

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

void EQsweep(qaqmc_local &qmc,int mstep){
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
	}
}

void MCsweep(qaqmc_local &qmc,int mstep,double &ma,double &m2,double &m4){
	ma = m2 = m4 = 0;
	int N = qmc.get_N();
	int Nb = qmc.get_Nb();
	std::vector<signed char> spins(N,0);
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		qmc.midpoint(spins.begin());
		double mtemp = 0;
		for(int i=0;i<N;i++){
			mtemp+=spins[i];
		}

		ma += std::abs(mtemp);
		m2 += std::pow(mtemp,2);
		m4 += std::pow(mtemp,4);
	}

	ma /= (mstep+1.1e-15);
	m2 /= (mstep+1.1e-15);
	m4 /= (mstep+1.1e-15);
}


double sfunc(double t){
	// std::cout << t << std::endl;
	return t;
}

int main(int argc, char const *argv[])
{
	int L = 10;
	int d = 1;
	int M = 10;
	int N = int(std::pow(double(L),d));

	std::vector<signed char> bc;

	for(int i=0;i<N;i++){
		bc.push_back(1);
	}

	std::vector<int> bst;

	sq_lattice(L,d,bst);

	double eJ,ma,m2,m4;

	// proj_local qmc(2*M,N,d*N,bst.data(),S,1,1,bc,bc);
	qaqmc_local qmc(2*M,N,d*N,bst.data(),&sfunc,1,1,bc,bc);
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


