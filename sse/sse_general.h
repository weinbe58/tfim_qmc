#ifndef __sse_general_INCLUDED__
#define __sse_general_INCLUDED__

#include <iostream>
#include "base/base.h"
#include "sse/sse_base.h"

// updates for non-uniform coupling strength on ising and tranverse field. 

class sse_general : public sse_base{
	protected:
		const int * bst;
		const double * Jb;
		const int Nb;
		const double Wh;
		double W,Wb;

		double beta;
		std::vector<double> P_op;
		void add_op(int);
		void remove_op(int);

	public:
		sse_general(double,const int, const int, const int[],const double[],const double);
		~sse_general() {};
		int get_Nb() {return Nb;}

};

sse_general::sse_general(double _beta,const int _N,const int _Nb, const int _bst[], const double _Jb[], const double _S) 
: sse_base::sse_base(_N), Nb(_Nb), bst(_bst), Jb(_Jb), beta(_beta), Wh(_N*(1-_S)), {
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
	Wb = 0.0;
	for(int i=0;i<Nb;i++){
		Wb += std::abs(_Jb[i]);
		P_op.push_back(Wb);
	}
	for(int i=0;i<Nb;i++){
		P_op[i] /= Wb;
	}

	Wb *= 2;
	W = _S*Wb + Wh;
}


void sse_general::add_op(int p){
	if(base::ran()*(base::M - sse_base::Nop + beta*W) < beta*W){
		if(base::ran()*W<Wh){
			base::opstr[p].o1=-1; 
			base::opstr[p].o2=std::floor(base::ran()*base::N);
			sse_base::Nop++;
		}
		else{
			auto iter = std::upper_bound(P_op.begin(), P_op.end(), ran());
			int ib = int(iter - P_op.begin()) << 1;
			int i = bst[ib];
			int j = bst[ib+1];
			if(Jb[ib>>1]*base::sP[i]*base::sP[j]<0){
				base::opstr[p].o1=i;
				base::opstr[p].o2=j;
				break;
			}
		}
	}
}

void sse_general::remove_op(int p){
	if(base::ran()*(base::M - sse_base::Nop + 1 + beta*W) < (base::M - sse_base::Nop + 1)){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
		sse_base::Nop--;
	}
}



#endif