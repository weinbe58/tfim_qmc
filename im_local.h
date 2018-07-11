#ifndef __im_local_INCLUDED__
#define __im_local_INCLUDED__

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include "base.h"
#include "im_base.h"



class im_local : public im_base{
	protected:
		const int * bst;
		const int Nb;

		void add_op_left(optype&,double);
		void add_op_right(optype&,double);
		void remove_op_left(optype&,double);
		void remove_op_right(optype&,double);
		double weight(optype,double);

	public:
		im_local(const int, const int, const int[], double,
		const double *,const int *, anneal_func, const int,const int,
		const std::vector<signed char>,const std::vector<signed char>);
		im_local(const int, const int, const int[], double,
		const double *,const int *, anneal_func, 
		const std::vector<signed char>,const std::vector<signed char>,
		const std::vector<signed char>,const std::vector<signed char>);
		~im_local() {};
		int get_Nb() {return Nb;}

};


im_local::im_local(const int _N,const int _Nb,const int _bst[],double _t_f,
	const double *_rpar,const int *_ipar,anneal_func func, const int _Fl,const int _Fr,
	const std::vector<signed char> _sL,const std::vector<signed char> _sR)
	 : im_base::im_base(_N,_t_f,_rpar,_ipar,func,_Fl,_Fr,_sL,_sR), Nb(_Nb), bst(_bst)
{
	for(int i=0;i<2*Nb;i++){
		if(bst[i]<0 || bst[i]>=N){
			std::cout << "bond index out of bounds" << std::endl;
			exit(-2);
		}
	}
}

im_local::im_local(const int _N,const int _Nb,const int _bst[],double _t_f,
	const double *_rpar,const int *_ipar,anneal_func func, 
	const std::vector<signed char> _Fl,const std::vector<signed char> _Fr,
	const std::vector<signed char> _sL,const std::vector<signed char> _sR)
	 : im_base::im_base(_N,_t_f,_rpar,_ipar,func,_Fl,_Fr,_sL,_sR), Nb(_Nb), bst(_bst)
{
	for(int i=0;i<2*Nb;i++){
		if(bst[i]<0 || bst[i]>=N){
			std::cout << "bond index out of bounds" << std::endl;
			exit(-2);
		}
	}
}



double im_local::weight(optype op,double t){
	if(op.o2>0){
		double S = S_func(t,im_base::t_f,im_base::rpar,im_base::ipar);
		if(op.o1>=0){
			return 2*S;
		}
		else{
			return ((1-S)*base::N)/Nb;
		}
	}
	return 1;
}

void im_local::add_op_left(optype &op,double t){
	double S = S_func(t,im_base::t_f,rpar,ipar);
	double W = 2*Nb*S + base::N*(1-S);
	
	if(base::ran()*(im_base::M_left - im_base::Nop_left + im_base::t_f*W) < im_base::t_f*W){
		double Wh = base::N*(1-S);
		if(base::ran()*W<Wh){
			op.o1=-1; 
			op.o2=std::floor(base::ran()*base::N);
			im_base::Nop_left++;
		}
		else{
			int ib = 2*std::floor(base::ran()*Nb);
			int i = bst[ib];
			int j = bst[ib+1];
			if(base::sP[i]*base::sP[j]>0){
				op.o1=i;
				op.o2=j;
				im_base::Nop_left++;
			}
		}
	}
}



void im_local::add_op_right(optype &op,double t){
	double S = S_func(t,im_base::t_f,rpar,ipar);
	double W = 2*Nb*S + base::N*(1-S);

	if(base::ran()*(im_base::M_right - im_base::Nop_right + im_base::t_f*W) < im_base::t_f*W){
		double Wh = base::N*(1-S);
		if(base::ran()*W<Wh){
			op.o1 = -1; 
			op.o2 = std::floor(base::ran()*base::N);
			im_base::Nop_right++;
		}
		else{
			int ib = 2*std::floor(base::ran()*Nb);
			int i = bst[ib];
			int j = bst[ib+1];
			if(base::sP[i]*base::sP[j]>0){
				op.o1=i;
				op.o2=j;
				im_base::Nop_right++;
			}
		}
	}
}

void im_local::remove_op_left(optype &op,double t){
	double S = S_func(t,im_base::t_f,rpar,ipar);
	double W = 2*Nb*S + base::N*(1-S);
	if(base::ran()*(im_base::M_left - im_base::Nop_left + 1 + im_base::t_f*W) < (im_base::M_left - im_base::Nop_left + 1)){
		op.o1=-1;
		op.o2=-1;
		im_base::Nop_left--;
	}
}

void im_local::remove_op_right(optype &op,double t){
	double S = S_func(t,im_base::t_f,rpar,ipar);
	double W = 2*Nb*S + base::N*(1-S);
	if(base::ran()*(im_base::M_right - im_base::Nop_right + 1 + im_base::t_f*W) < (im_base::M_right - im_base::Nop_right + 1)){
		op.o1=-1;
		op.o2=-1;
		im_base::Nop_right--;
	}
}




#endif