#ifndef __sse_base_INCLUDED__
#define __sse_base_INCLUDED__

#include "base.h"

class sse_base : base
{
	protected:
		int Nop;

		virtual bool add_op(int) = 0;
		virtual bool remove_op(int) = 0;

	public:
		sse_base(const int _N) : base::base::(_N,_N), Nop(0) {}
		~sse_base();

		void check_M(void);
	
};

void sse_base::diagonal_update(){
	for(i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];
	}// end for(int i=1;...

	for(int p=0;p<base::M;p++){
		if(base::opstr[p].o2<0){
			add_op(p);
		}
		else if(opstr_l[p].o1>-2){// try to remove diagonal operator
			remove_op(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}// end (opstr_l[p].o2<0)
	}
}

void sse_base::check_M(){
	int Mnew=(Nop*4)/3;
	if(Mnew > base::M){	
		std::vector<optype> opstr_temp(base::opstr);
		base::opstr.resize(Mnew);
		for(int j=0;j<base::M;j++){ base::opstr[j]=base::opstr[j];}
		for(int j=base::M;j<Mnew;j++){base::opstr[j].o1=-1; base::opstr[j].o2=-1;}
		base::M=Mnew;
	}
}

#endif