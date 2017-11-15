#ifndef __sse_local_INCLUDED__
#define __sse_local_INCLUDED__

#include <iostream>
#include "base.h"
#include "sse_base.h"

class sse_local : public sse_base{
	protected:
		const int * bst;
		const int Nb;
		const double W,WJ,Wh;
		double beta;

		void add_op(int);
		void remove_op(int);

	public:
		sse_local(const int, const int, const int[],const double, const double);
		~sse_local() {};

};

sse_local::sse_local(const int _N,const int _Nb, const int _bst[], const double _S, const double _beta) 
: sse_base::sse_base(_N), Nb(_Nb), bst(_bst), beta(_beta), WJ(2*_Nb*(1-_S)), Wh(_N*_S), W(2*_Nb*(1-_S)+_N*_S) {
	for(int i=0;i<2*Nb;i++){
		if(bst[i]<0 || bst[i]>=N){
			std::cout << "bond index out of bounds" << std::endl;
			exit(-2);
		}
	}
	if(_S<0.0 || _S>1.0){
		std::cout << "S must be in [0,1]." << std::endl;
		exit(-3);
	}
}


void sse_local::add_op(int p){
	if(base::ran()*(base::M - sse_base::Nop + beta*W) < beta*W){
		if(base::ran()*(Wh+WJ)<Wh){
			base::opstr[p].o1=-1; 
			base::opstr[p].o2=std::floor(base::ran()*base::N);
			sse_base::Nop++;
		}
		else{
			int ib = 2*std::floor(base::ran()*Nb);
			int i = bst[ib];
			int j = bst[ib+1];
			if(base::sP[i]*base::sP[j]>0){
				base::opstr[p].o1=i;
				base::opstr[p].o2=j;
				sse_base::Nop++;
			}
		}
	}
}

void sse_local::remove_op(int p){
	if(base::ran()*(base::M - sse_base::Nop + 1 + beta*W) < (base::M - sse_base::Nop + 1)){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
		sse_base::Nop--;
	}
}



#endif