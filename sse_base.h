#ifndef __sse_base_INCLUDED__
#define __sse_base_INCLUDED__

#include "base.h"

class sse_base : public base{
	protected:
		int Nop;
		virtual void add_op(int) = 0;
		virtual void remove_op(int) = 0;

	public:
		sse_base(const int _N);
		~sse_base() {};

		void diagonal_update(void);
		void check_M(void);
};


sse_base::sse_base(const int _N) : base::base(_N,_N), Nop(0) {
	for(int p=0;p<base::M;p++){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
	}
}

void sse_base::diagonal_update(){
	// this->print_opstr(false);

	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];
		// std::cout << int(base::sL[i]) << "  ";
	}// end for(int i=1;...
	// std::cout << std::endl;

	for(int p=0;p<base::M;p++){
		if(base::opstr[p].o2<0){
			add_op(p);
		}
		else if(base::opstr[p].o1>-2){// try to remove diagonal operator
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
		base::X.resize(4*Mnew);
		for(int j=0;j<base::M;j++){ base::opstr[j]=base::opstr[j];}
		for(int j=base::M;j<Mnew;j++){base::opstr[j].o1=-1; base::opstr[j].o2=-1;}
		base::M=Mnew;
	}
}

#endif