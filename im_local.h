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

		void add_op_left(int);
		void add_op_right(int);
		void remove_op_left(int);
		void remove_op_right(int);

	public:
		im_local(const int, const int, const int[], double,
		const double *,const int *, anneal_func, const int,const int,
		const std::vector<signed char>,const std::vector<signed char>);
		// im_local(const int, const int, const int[], const double,
		// const double *,const int *, anneal_func, const int,const int);
		~im_local() {};
		int get_Nb() {return Nb;}
		int update_times(int);

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


void im_local::add_op_left(int p){
	double S = S_func(im_base::times[p],rpar,ipar);
	double W = 2*Nb*S + base::N*(1-S);
	
	if(base::ran()*(im_base::M_left - im_base::Nop_left + im_base::t_f*W) < im_base::t_f*W){
		double Wh = base::N*(1-S);
		if(base::ran()*W<Wh){
			base::opstr[p].o1=-1; 
			base::opstr[p].o2=std::floor(base::ran()*base::N);
			im_base::Nop_left++;
		}
		else{
			int ib = 2*std::floor(base::ran()*Nb);
			int i = bst[ib];
			int j = bst[ib+1];
			if(base::sP[i]*base::sP[j]>0){
				base::opstr[p].o1=i;
				base::opstr[p].o2=j;
				im_base::Nop_left++;
			}
		}
	}
}



void im_local::add_op_right(int p){
	double S = S_func(im_base::times[p],rpar,ipar);
	double W = 2*Nb*S + base::N*(1-S);

	if(base::ran()*(im_base::M_right - im_base::Nop_right + im_base::t_f*W) < im_base::t_f*W){
		double Wh = base::N*(1-S);
		if(base::ran()*W<Wh){
			base::opstr[p].o1=-1; 
			base::opstr[p].o2=std::floor(base::ran()*base::N);
			im_base::Nop_right++;
		}
		else{
			int ib = 2*std::floor(base::ran()*Nb);
			int i = bst[ib];
			int j = bst[ib+1];
			if(base::sP[i]*base::sP[j]>0){
				base::opstr[p].o1=i;
				base::opstr[p].o2=j;
				im_base::Nop_right++;
			}
		}
	}
}

void im_local::remove_op_left(int p){
	double S = S_func(im_base::times[p],rpar,ipar);
	double W = 2*Nb*S + base::N*(1-S);
	if(base::ran()*(im_base::M_left - im_base::Nop_left + 1 + im_base::t_f*W) < (im_base::M_left - im_base::Nop_left + 1)){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
		im_base::Nop_left--;
	}
}

void im_local::remove_op_right(int p){
	double S = S_func(im_base::times[p],rpar,ipar);
	double W = 2*Nb*S + base::N*(1-S);
	if(base::ran()*(im_base::M_right - im_base::Nop_right + 1 + im_base::t_f*W) < (im_base::M_right - im_base::Nop_right + 1)){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
		im_base::Nop_right--;
	}
}

int im_local::update_times(int n){
	// for(auto t: times){
	// 	std::cout << std::setw(15) << t;
	// }
	// std::cout << std::endl;
	
	n = std::min(n,std::min(im_base::M_left,im_base::M_right));
	n = std::max(1,n);
	int update=0;
	int m = std::floor(base::ran()*(im_base::M_left-n));

	double t_1 = im_base::times[m];
	double t_2 = im_base::times[m+n];
	std::vector<double> new_times;
	for(int p=0;p<n;p++)
		new_times.push_back(t_1+base::ran()*(t_2-t_1));

	std::sort(new_times.begin(),new_times.end(),std::greater<double>());

	double W_new=1;
	double W_old=1;

	for(int p=0;p<n;p++){
		double S_new = S_func(new_times[p],rpar,ipar);
		double S_old = S_func(im_base::times[m+p],rpar,ipar);
		if(base::opstr[p].o2>0){
			if(base::opstr[p].o1>0){
				W_new *= 2*Nb*S_new;
				W_old *= 2*Nb*S_old;
			}
			else{
				W_new *= base::N*S_new;
				W_old *= base::N*S_old;				
			}
		}
	}

	if(base::ran()*W_old<W_new){
		std::copy(new_times.begin(),new_times.end(),im_base::times.begin()+m);
		update++;
	}


	m = M_left+std::floor(base::ran()*(im_base::M_right-n));

	t_1 = im_base::im_base::times[m];
	t_2 = im_base::im_base::times[m+n];

	for(int p=0;p<n;p++)
		new_times[p]=(t_1+base::ran()*(t_2-t_1));

	std::sort(new_times.begin(),new_times.end(),std::less<double>());

	W_new=1;
	W_old=1;

	for(int p=0;p<n;p++){
		double S_new = S_func(new_times[p],rpar,ipar);
		double S_old = S_func(im_base::im_base::times[m+p],rpar,ipar);
		if(base::opstr[p].o2>0){
			if(base::opstr[p].o1>0){
				W_new *= 2*Nb*S_new;
				W_old *= 2*Nb*S_old;
			}
			else{
				W_new *= base::N*S_new;
				W_old *= base::N*S_old;				
			}
		}
	}
	// std::cout << W_old << std::setw(20) << W_new << std::endl;

	if(base::ran()*W_old<W_new){
		std::copy(new_times.begin(),new_times.end(),im_base::times.begin()+m);
		update++;
	}

	// for(auto t: times){
	// 	std::cout << std::setw(15) << t;
	// }
	// std::cout << std::endl;

	return update;
}



#endif