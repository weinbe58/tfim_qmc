#ifndef __im_local_INCLUDED__
#define __im_local_INCLUDED__

#include <iostream>
#include "base.h"
#include "im_base.h"

class im_local : public im_base{
	protected:
		const int * bst;
		const int Nb;
		const double W,Wh;
		double beta;

		void add_op(int);
		void remove_op(int);

	public:
		im_local(const int, const int, const int[], const double, const int, const int,
		const std::vector<signed char>,const std::vector<signed char>);

		~im_local() {};
		int get_Nb() {return Nb;}

};

im_local::im_local(const int _N,const int _Nb, const int _bst[], const double _S) 
: im_base::im_base(_N), Nb(_Nb), bst(_bst), beta(_beta), Wh(_N*(1-_S)), W(2*_Nb*_S+_N*(1-_S)) {
	for(int i=0;i<2*Nb;i++){
		if(bst[i]<0 || bst[i]>=_N){
			std::cout << "bond index out of bounds" << std::endl;
			exit(-2);
		}
	}
	if(_S<0.0 || _S>1.0){
		std::cout << "S must be in [0,1]." << std::endl;
		exit(-3);
	}
}


void im_local::add_op(int p){
	if(base::ran()*(base::M - im_base::Nop + beta*W) < beta*W){
		if(base::ran()*W<Wh){
			base::opstr[p].o1=-1; 
			base::opstr[p].o2=std::floor(base::ran()*base::N);
			im_base::Nop++;
		}
		else{
			int ib = 2*std::floor(base::ran()*Nb);
			int i = bst[ib];
			int j = bst[ib+1];
			if(base::sP[i]*base::sP[j]>0){
				base::opstr[p].o1=i;
				base::opstr[p].o2=j;
				im_base::Nop++;
			}
		}
	}
}

void im_local::remove_op(int p){
	if(base::ran()*(base::M - im_base::Nop + 1 + beta*W) < (base::M - im_base::Nop + 1)){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
		im_base::Nop--;
	}
}



#endif