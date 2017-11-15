#ifndef __sse_base_INCLUDED__
#define __sse_base_INCLUDED__

#include "base.h"


class sse_base : base
{
	protected:
		const int * bst;
		const int Nb;
		int Nop;


	virtual bool add_op(int) = 0;
	virtual bool remove_op(int) = 0;
	virtual int weight1(void) = 0;
	virtual int weight2(int) = 0;

	public:
		sse_base(const int, const int, const * int);
		~sse_base();

		void check_M(void);
	
};


sse_base::sse_base(const int _N,
		 const int _Nb,
		 const int * _bst,
		 const std::vector<signed char> _sL) : base::base::(_N,_N), Nb(_Nb), bst(_bst) {

}

void sse_base::diagonal_update(){
	for(i=0;i<N;i++){
		spins[i]=spinsL[i];
	}// end for(int i=1;...

	for(int p=0;p<M;p++){
		if(opstr[p].o2<0){
			if(add_op(p)){
				int i = weight1();
				int j = weight2(i);
				if(i == j){ 
					opstr[p].o1=-1; 
					opstr[p].o2=i;
					Nop++;
				}
				else if(spins[i]*spins[j]*J[i][j]<0){
					opstr[p].o1=i;
					opstr[p].o2=j;
					Nop++;
				}
			}
		}
		else if(opstr_l[p].o1>-2){// try to remove diagonal operator
			if(remove_op(p)){
				opstr[p].o1=-1;
				opstr[p].o2=-1;
				Nop--;
			}
		}
		else{
			spins[ opstr[p].o2 ] *= -1;
		}// end (opstr_l[p].o2<0)
	}
}

void sse_base::check_M(){
	int Mnew=(Nop*4)/3;
	if(Mnew > M){	
		std::vector<optype> opstr_temp(opstr);
		opstr.resize(Mnew);
		for(int j=0;j<M;j++){ opstr[j]=opstr[j];}
		for(int j=M;j<Mnew;j++){opstr[j].o1=-1; opstr[j].o2=-1;}
		M=Mnew;
	}
}

#endif