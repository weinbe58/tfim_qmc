#ifndef __qaqmc_local_INCLUDED__
#define __qaqmc_local_INCLUDED__

#include <iostream>
#include <algorithm>
#include <cmath>
#include "base.h"
#include "proj_base.h"

class qaqmc_glass : public proj_base{
	protected:
		const int * bst;
		const double * Jb;
		int * ipar;
		double * rpar;
		const int Nb;
		const double Smax;
		double beta;
		std::vector<double> P_op;
		double Wb;
		double (*S_func)(double,double*,int*);
		void move_op(int);

	public:
		qaqmc_glass(int,const int,const int, const int[],const double[],double*,int*,double (*)(double,double*,int*), const int, const int,
			const short_vec,const short_vec);
		qaqmc_glass(int,const int,const int, const int[],const double[],double*,int*,double (*)(double,double*,int*), const int, const int);

		~qaqmc_glass() {};
		int get_Nb() {return Nb;}
		

};


qaqmc_glass::qaqmc_glass(int _M,const int _N,const int _Nb, const int _bst[],const double _Jb[],double *_rpar,int *_ipar,double (* func)(double,double*,int*),
		 	const int _Fl, const int _Fr, const short_vec _sL,const short_vec _sR) :
proj_base::proj_base(_M,_N,_Fl,_Fr,_sL,_sR), rpar(_rpar), ipar(_ipar), Nb(_Nb), bst(_bst), Jb(_Jb), S_func(func), Smax(func(1.0,_rpar,_ipar)){
	for(int i=0;i<2*Nb;i++){
		if(bst[i]<0 || bst[i]>=_N){
			std::cout << "bond index out of bounds" << std::endl;
			exit(-2);
		}
	}

	Wb = 0.0;
	for(int i=0;i<Nb;i++){
		Wb += std::abs(Jb[i]);
		P_op.push_back(Wb);
	}
	for(int i=0;i<Nb;i++){
		P_op[i] /= Wb;
	}

	Wb *= 2;
	
	
}

qaqmc_glass::qaqmc_glass(int _M,const int _N,const int _Nb, const int _bst[],const double _Jb[],double *_rpar,int *_ipar,double (* func)(double,double*,int*),
		 	const int _Fl, const int _Fr) : proj_base::proj_base(_M,_N,_Fl,_Fr),
		 	rpar(_rpar), ipar(_ipar), Nb(_Nb), bst(_bst), Jb(_Jb), S_func(func), Smax(func(1.0,_rpar,_ipar)){

	for(int i=0;i<2*Nb;i++){
		if(bst[i]<0 || bst[i]>=_N){
			std::cout << "bond index out of bounds" << std::endl;
			exit(-2);
		}
	}

	Wb = 0.0;
	for(int i=0;i<Nb;i++){
		Wb += std::abs(Jb[i]);
		P_op.push_back(Wb);
	}
	for(int i=0;i<Nb;i++){
		P_op[i] /= Wb;
	}

	Wb *= 2;
	
	
}


void qaqmc_glass::move_op(int p){
	double S;

	if(p < base::M/2){
		S = S_func(2*double(p+1)/base::M,rpar,ipar);
	}
	else{
		S = Smax-S_func(2*double(p-base::M/2)/base::M,rpar,ipar);
	}
	double W = Wb*S + base::N*(1-S);
	while(true){
		if(base::ran()*W < Wb*S ){
			auto iter = std::upper_bound(P_op.begin(), P_op.end(), ran());
			int ib = int(iter - P_op.begin()) << 1;
			int i = bst[ib];
			int j = bst[ib+1];
			// std::cout << ib << "  " << int(base::sP[i]) << "  " << int(base::sP[j]) << std::endl; 
			if(Jb[ib>>1]*base::sP[i]*base::sP[j]<0){
				base::opstr[p].o1=i;
				base::opstr[p].o2=j;
				break;
			}
		}
		else{
			base::opstr[p].o1 = -1; 
			base::opstr[p].o2 = std::floor(base::ran()*base::N);
			break;

		}
	}
}




#endif
