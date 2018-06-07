

#include "im_local.h"

#include <iostream>
#include <iomanip>
#include <utility>
#include <string>



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

int EQsweep(im_local &qmc,int mstep){
	int n_up = 0;
	int n = qmc.get_M();
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		qmc.check_M();
		// std::cout << qmc.get_M() << std::endl;
		n_up += qmc.update_times(n);
		std::cout << n << std::setw(10) << n_up << std::endl;
		if(i%50==0){
			if(n_up<50){
				if(n < qmc.get_M())
					n++;
			}
			else{
				if(n > 1)
					n--;
			}
			n_up = 0;
		}
	}

}

std::pair<double,double> MCsweep(im_local &qmc,int mstep,int n_up){
	double m2=0;
	double m4=0;
	std::vector<signed char> spins(qmc.get_N(),0);
	for(int i=0;i<mstep;i++){
		qmc.diagonal_update();
		qmc.cluster_update();
		qmc.update_times(n_up);

		qmc.midpoint(spins.begin());
		double mtemp = 0;
		for(auto s:spins){mtemp += s;}
		m2 += std::pow(mtemp,2);
		m4 += std::pow(mtemp,4);
	}

	m2 /= (mstep+1.1e-15);
	m4 /= (mstep+1.1e-15);

	return std::make_pair(m2,m4);
}


double sfunc(double t,const double * rpar,const int * ipar){
	return t;
}

int main(int argc, char const *argv[])
{
	int L=4;
	int mstep=10000;
	int nbin=1000;
	double t_f=1.0;
	std::string outfile;
	// std::cin >> L >> t_f >> nbin >> mstep >> outfile;

	int d = 1;
	int N = std::pow(L,d);


	std::vector<signed char> bc;

	for(int i=0;i<N;i++){
		bc.push_back(1);
	}

	std::vector<int> bst;

	sq_lattice(L,d,bst);

	// std::fstream fs(outfile.c_str(), std::fstream::app);

	// if(!fs.is_open()){
	// 	std::cout << "failed to open file!" << std::endl;
	// 	std::exit(-1);
	// }
	// std::ostream * buffer = &fs;
	std::ostream * buffer = &std::cout;

	double S_f = 1.0;
	im_local qmc(N,d*N,bst.data(),t_f,&S_f,NULL,&sfunc,1,1,bc,bc);
	int n_up = EQsweep(qmc,mstep);

	// for(int i=0;i<nbin;i++){
	// 	std::cout << "bin: " << i << std::endl;
	// 	std::pair<double,double> pair = MCsweep(qmc,mstep,n_up);
	// 	(*buffer) << std::scientific << std::setprecision(10) << std::setw(20);
	// 	(*buffer) << std::get<0>(pair)/(N*N) << std::setw(20);
	// 	(*buffer) << std::get<1>(pair)/(N*N*N*N) << std::setw(20);
	// 	(*buffer) << std::endl;
	// }


	return 0;
}

