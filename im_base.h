#ifndef __im_base_INCLUDED__
#define __im_base_INCLUDED__

#include "base.h"

class im_base : public base{
	protected:
		int Nop_left,Nop_right;
		int M_left,M_right;
		virtual void add_op(int) = 0;
		virtual void remove_op(int) = 0;

	public:
		im_base(const int _N);
		~im_base() {};

		void diagonal_update(void);
		void check_M(void);
};


im_base::im_base(const int _N,
				 const int _Fl,
				 const int _Fr,
				 const std::vector<signed char> _sL,
				 const std::vector<signed char> _sR) : base::base(2*_N,_N,_Fl,_Fr,_sL,_sR),
				 Nop_left(0), Nop_right(0), M_left(_N), M_right(_N)
{
	for(int p=0;p<base::M;p++){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=std::floor(base::N*base::ran());
	}
}


void im_base::diagonal_update(){
	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];		
	}

	for(int p=0;p<M_left;p++){
		if(base::opstr[p].o2 == -1){
			add_op(p);
		}
		else if(base::opstr[p].o1>-2){
			remove_op(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}
	}

	for(int p=M_left;p<M_right;p++){
		if(base::opstr[p].o2 == -1){
			add_op(p);
		}
		else if(base::opstr[p].o1>-2){
			remove_op(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}
	}
}

void im_base::check_M(){
	int M_left_new=(Nop_left*4)/3;
	if(M_left_new > M_left){
		int nr = M_left_new - M_left;
		optype val;
		base::opstr.insert(base::opstr.begin(),val,nl);
	}


	int M_right_new=(Nop_right*4)/3;
	if(M_right_new > M_right){	
		int nl = M_right_new - M_right;
		optype val;
		base::opstr.insert(base::opstr.end(),val,nr);
	}
	base:M = base::opstr.size();
}

#endif