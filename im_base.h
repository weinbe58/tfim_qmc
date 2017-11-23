#ifndef __im_base_INCLUDED__
#define __im_base_INCLUDED__

#include "base.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>


typedef double (* anneal_func)(double,const double*,const int*);

class im_base : public base{
	protected:
		int Nop_left,Nop_right;
		int M_left,M_right;
		std::vector<double> times;
		double t_f;
		const int * ipar;
		const double * rpar;

		anneal_func S_func;

		virtual void add_op_left(int) = 0;
		virtual void add_op_right(int) = 0;
		virtual void remove_op_left(int) = 0;
		virtual void remove_op_right(int) = 0;

	public:
		im_base(const int, const double,const double *,const int *, anneal_func,
		const int,const int,const std::vector<signed char>,const std::vector<signed char>);
		// im_base(const int, const double,const double *,const int *, anneal_func,
		// const int,const int);

		~im_base() {};

		void diagonal_update(void);
		virtual int update_times(int)=0;
		void check_M(void);
		void midpoint(std::vector<signed char>::iterator);
};


im_base::im_base(const int _N,
				 double _t_f,
				 const double *_rpar,
				 const int *_ipar,
				 anneal_func func,
				 const int _Fl,
				 const int _Fr,
				 const std::vector<signed char> _sL,
				 const std::vector<signed char> _sR) : base::base(2*_N,_N,_Fl,_Fr,_sL,_sR),
				 S_func(func),rpar(_rpar),ipar(_ipar),t_f(_t_f),Nop_left(0), Nop_right(0), M_left(_N), M_right(_N)
{
	for(int i=0;i<2*_N;i++)
		times.push_back(base::ran()*_t_f);

	std::sort(times.begin(),times.begin()+_N,std::greater<double>());
	std::sort(times.begin()+_N,times.end(),std::less<double>());

	for(int p=0;p<base::M;p++){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
	}
}


void im_base::diagonal_update(){
	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];		
	}
	
	for(int p=0;p<M_left;p++){
		
		if(base::opstr[p].o2 == -1){
			add_op_left(p);
		}
		else if(base::opstr[p].o1>-2){
			remove_op_left(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}
	}

	for(int p=M_left;p<base::M;p++){
		if(base::opstr[p].o2 == -1){
			add_op_right(p);
		}
		else if(base::opstr[p].o1>-2){
			remove_op_right(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}
	}
}





void im_base::midpoint(std::vector<signed char>::iterator spins){
	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];
	}

	int p = 0;
	auto op = base::opstr.begin();

	while(p<=M_left){
		if(op->o1==-2){
			base::sP[ op->o2 ] *= -1;
		}
		op++;
		p++;
	}

	std::copy(base::sP.begin(),base::sP.end(),spins);
	
}


void im_base::check_M(){
	int M_left_new=(Nop_left*4)/3;
	
	if(M_left_new > M_left){
		int nl = M_left_new - M_left;
		optype val;
		base::opstr.insert(base::opstr.begin(),nl,val);
		double t_min = times[0];
		for(int i=0;i<nl;i++){
			double t = base::ran()*t_min;
			times.insert(times.end(),t);
		}
		std::sort(times.begin(),times.begin()+nl,std::greater<double>());
		M_left = M_left_new;
	}

	int M_right_new=(Nop_right*4)/3;
	if(M_right_new >= M_right){	
		int nr = M_right_new - M_right;
		optype val;
		base::opstr.insert(base::opstr.end(),nr,val);
		double t_min = times[M_right-1];
		for(int i=0;i<nr;i++){
			double t = base::ran()*t_min;
			times.insert(times.end(),t);
		}
		std::sort(times.end()-nr,times.end(),std::less<double>());
		M_right = M_right_new;
	}
	base:M = base::opstr.size();
	X.resize(4*base::M);
	// std::cout << base::opstr.size() << std::setw(10) << times.size() << std::endl;
}

#endif