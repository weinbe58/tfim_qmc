#ifndef __qaqmc_local_INCLUDED__
#define __qaqmc_local_INCLUDED__

#include <iostream>
#include "base.h"
#include "proj_base.h"

class qaqmc_local : public proj_base{
	protected:
		const int * bst;
		const int Nb;
		const double Smax;
		double beta;
		double (*S_func)(double);
		void move_op(int);

	public:
		qaqmc_local(int,const int,const int, const int[], double (*)(double), const int, const int,
			const std::vector<signed char>,const std::vector<signed char>);
		~qaqmc_local() {};
		int get_Nb() {return Nb;}

};


qaqmc_local::qaqmc_local(int _M,const int _N,const int _Nb, const int _bst[],double (* func)(double),
		 const int _Fl, const int _Fr, const std::vector<signed char> _sL,const std::vector<signed char> _sR) :
proj_base::proj_base(_M,_N,_Fl,_Fr,_sL,_sR), Nb(_Nb), bst(_bst), S_func(func), Smax(func(1.0)){
	for(int i=0;i<2*Nb;i++){
		if(bst[i]<0 || bst[i]>=_N){
			std::cout << "bond index out of bounds" << std::endl;
			exit(-2);
		}
	}
}


void qaqmc_local::move_op(int p){
	double S;
	if(p < base::M/2){
		S = S_func(2*double(p+1)/base::M);
	}
	else{
		S = Smax-S_func(2*double(p-base::M/2)/base::M);
	}
	double W = 2*Nb*S + base::N*(1-S);;
	while(true){
		if(base::ran()*W< base::N*(1-S) ){
			base::opstr[p].o1=-1; 
			base::opstr[p].o2=std::floor(base::ran()*base::N);
			break;
		}
		else{
			int ib = 2*std::floor(base::ran()*Nb);
			int i = bst[ib];
			int j = bst[ib+1];
			if(base::sP[i]*base::sP[j]>0){
				base::opstr[p].o1=i;
				base::opstr[p].o2=j;
				break;
			}
		}
	}
}




#endif