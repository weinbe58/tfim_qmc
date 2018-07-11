

#include "im_local.h"

#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <functional>



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

template<class qmc_class>
void EQsweep_1(qmc_class &qmc,int mstep){
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		qmc.check_M();
	}
}

template<class qmc_class>
int EQsweep_2(qmc_class &qmc,int mstep){
	int n_up = 0;
	int n = 1;
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		double r_accept = qmc.update_times(n,10);

		if(r_accept>0.5){
			if(n < qmc.get_M())
				n++;
		}
		else{
			if(n > 1)
				n--;
		}
		n_up += n;
	}

	return n_up/mstep;
}

template<class qmc_class>
std::pair<double,double> MCsweep(qmc_class &qmc,int mstep,int n_up){
	double m2=0;
	double m4=0;
	std::vector<signed char> spins(qmc.get_N(),0);
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		qmc.update_times(n_up,10);
		
		qmc.midpoint(spins.begin());
		double mtemp = std::accumulate(spins.begin(),spins.end(),0);
		
		m2 += std::pow(mtemp,2);
		m4 += std::pow(mtemp,4);
	}
	// qmc.print_opstr(false);
	m2 /= (mstep+1.1e-15);
	m4 /= (mstep+1.1e-15);

	return std::make_pair(m2,m4);
}


double sfunc(double t,double t_f,const double * rpar,const int * ipar){
	return rpar[0]*(t/t_f);
}

int main(int argc, char const *argv[])
{
	int L=10;
	int mstep=10000;
	int nbin=100;
	double t_f=0.5;
	double S_f = 0.5;
	std::string outfile;
	std::cin >> L >> t_f >> nbin >> mstep >> outfile;

	int d = 1;
	int N = std::pow(L,d);


	std::vector<signed char> bc;

	for(int i=0;i<N;i++){
		bc.push_back(1);
	}

	std::vector<int> bst;

	sq_lattice(L,d,bst);

	std::fstream fs(outfile.c_str(), std::fstream::app);

	if(!fs.is_open()){
		std::cout << "failed to open file!" << std::endl;
		std::exit(-1);
	}
	std::ostream * buffer = &fs;
	// std::ostream * buffer = &std::cout;

	im_local qmc(N,d*N,bst.data(),t_f,&S_f,NULL,&sfunc,1,1,bc,bc);
	
	EQsweep_1(qmc,mstep);
	int n_up = EQsweep_2(qmc,mstep);
	// int n_up = 1;

	for(int i=0;i<nbin;i++){
		std::cout << "bin: " << i+1 << std::endl;
		std::pair<double,double> pair = MCsweep(qmc,mstep,n_up);
		(*buffer) << std::scientific << std::setprecision(10);
		(*buffer) << std::setw(20) << t_f;
		(*buffer) << std::setw(20) << std::get<0>(pair)/(N*N);
		(*buffer) << std::setw(20) << std::get<1>(pair)/(N*N*N*N);
		(*buffer) << std::endl;
	}


	return 0;
}

