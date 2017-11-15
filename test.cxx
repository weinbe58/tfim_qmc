
#include "sse_local.h"
#include <iostream>
#include <iomanip>
#include <tuple>

void EQsweep(sse_local &qmc,int mstep){
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		qmc.check_M();
	}
}

void MCsweep(sse_local &qmc,int mstep,std::vector<int> bst,int &ef,int &et,
				double &eJ,double &ma,double &m2,double &m4){
	ef = et = 0;
	ma = m2 = m4 = eJ = 0;
	int M = qmc.get_M();
	int Nm = 0;
	int Nb = bst.size();
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		// qmc.print_opstr(false);
		for(auto tup=qmc.state_begin();tup != qmc.state_end();tup++){
			double msprop = 2*std::get<1>(*tup);
			double msprop2 = msprop*msprop;

			ma += std::abs(msprop);
			m2 += msprop2;
			m4 += msprop2*msprop2;
			for(int b=0;b<Nb;b+=2){
				int i = bst[b];
				int j = bst[b+1];
				eJ += std::get<2>(*tup)[i]*std::get<2>(*tup)[j];
				// std::cout << Nm <<"   "<< std::get<2>(*tup)[i]*std::get<2>(*tup)[j] << std::endl;
			}

			Nm++;
		}
	}

	ma /= (Nm+1.1e-15);
	m2 /= (Nm+1.1e-15);
	m4 /= (Nm+1.1e-15);
	eJ /= (Nm+1.1e-15);
}


int main(int argc, char const *argv[])
{
	const int L = 4;
	const int d = 1;
	const double S=0.3;
	const double beta = 10.0;
	int N = 1;

	for(int i=0;i<d;i++){
		N *= L;
	}


	std::vector<int> bst;

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

	double eJ,ma,m2,m4;
	int et,ef;

	sse_local qmc(N,d*N,bst.data(),S,beta);
	for(int i=0;i<5;i++){
		EQsweep(qmc,10000);
	}
	for(int i=0;i<5;i++){
		MCsweep(qmc,10000,bst,et,ef,eJ,ma,m2,m4);
		std::cout << std::scientific << std::setprecision(10) << std::setw(20);
		// std::cout << ma/N << std::setw(20);
		// std::cout << -N*S*double(ef)/(et+1.1e-15)-(1-S)*eJ << std::setw(20);
		std::cout << -N*S*double(ef)/(et+1.1e-15) << std::setw(20);
		std::cout << -(1-S)*eJ << std::setw(20);
		std::cout << m2/(N*N) << std::setw(20);
		std::cout << m4/(N*N*N*N) << std::setw(20);
		std::cout << std::endl;
	}


	return 0;
}


