#ifndef __proj_local_INCLUDED__
#define __proj_local_INCLUDED__

#include <iostream>
#include "base.h"
#include "proj_base.h"

class proj_local : public proj_base{
	protected:
		const int * bst;
		const int Nb;
		const double W,Wh;

		void move_op(int);

	public:
		proj_local(int, const int, const int, const int[],const double);
		proj_local(int,const int, const int, const int[], const double, const int, const int,
		const short_vec,const short_vec);
		~proj_local() {};

		int inline get_Nb() {return Nb;}	

};

proj_local::proj_local(int _M,const int _N,const int _Nb, const int _bst[], const double _S) 
: proj_base::proj_base(_M,_N), Nb(_Nb), bst(_bst), Wh(_N*(1-_S)), W(2*_Nb*_S+_N*(1-_S)) {
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


proj_local::proj_local(int _M,const int _N,const int _Nb, const int _bst[], const double _S, const int _Fl, const int _Fr,
		const short_vec _sL,const short_vec _sR) :
proj_base::proj_base(_M,_N,_Fl,_Fr,_sL,_sR), Nb(_Nb), bst(_bst), W(2*_Nb*_S+_N*(1-_S)), Wh(_N*(1-_S)){
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


void proj_local::move_op(int p){
	while(true){
		if(base::ran()*W<Wh){
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
