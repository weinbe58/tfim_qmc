#ifndef __im_base_INCLUDED__
#define __im_base_INCLUDED__

#include "base/base.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>


typedef double (* anneal_func)(double,double,const double*,const int*);

class im_base : public base{
	protected:
		int Nop_left,Nop_right;
		int M_left,M_right;
		std::vector<double> t_left,t_right;
		double t_f;
		const int * ipar;
		const double * rpar;

		anneal_func S_func;

		virtual void add_op_left(optype&,double) = 0;
		virtual void add_op_right(optype&,double) = 0;
		virtual void remove_op_left(optype&,double) = 0;
		virtual void remove_op_right(optype&,double) = 0;
		virtual double weight(optype,double) = 0;

		double inline time(int);

	public:
		im_base(const int, const double,const double *,const int *, anneal_func,
		const int,const int,const short_vec,const short_vec);
		im_base(const int, const double,const double *,const int *, anneal_func,
		const short_vec,const short_vec,
		const short_vec,const short_vec);


		~im_base() {};

		void diagonal_update(void);
		double update_times(int,int);
		void check_M(void);
		void midpoint(short_vec::iterator);
		int inline get_M_left(){return M_left;}
		int inline get_M_right(){return M_right;}
};


im_base::im_base(const int _N,
				 double _t_f,
				 const double *_rpar,
				 const int *_ipar,
				 anneal_func func,
				 const int _Fl,
				 const int _Fr,
				 const short_vec _sL,
				 const short_vec _sR) : base::base(2*_N,_N,_Fl,_Fr,_sL,_sR),
				 S_func(func),rpar(_rpar),ipar(_ipar),t_f(_t_f),Nop_left(0), Nop_right(0), M_left(_N), M_right(_N)
{
	for(int i=0;i<_N;i++)
		t_left.push_back(base::ran()*_t_f);

	for(int i=0;i<_N;i++)
		t_right.push_back(base::ran()*_t_f);

	std::sort(t_left.begin(),t_left.end());
	std::sort(t_right.begin(),t_right.end());

	for(int p=0;p<base::M;p++){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
	}
}

im_base::im_base(const int _N,
				 double _t_f,
				 const double *_rpar,
				 const int *_ipar,
				 anneal_func func,
				 const short_vec _Fr,
				 const short_vec _Fl,
				 const short_vec _sL,
				 const short_vec _sR) : base::base(2*_N,_N,_Fl,_Fr,_sL,_sR),
				 S_func(func),rpar(_rpar),ipar(_ipar),t_f(_t_f),Nop_left(0), Nop_right(0), M_left(_N), M_right(_N)
{
	for(int i=0;i<_N;i++)
		t_left.push_back(base::ran()*_t_f);

	for(int i=0;i<_N;i++)
		t_right.push_back(base::ran()*_t_f);

	std::sort(t_left.begin(),t_left.end());
	std::sort(t_right.begin(),t_right.end());

	for(int p=0;p<base::M;p++){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
	}
}


void im_base::diagonal_update(){
	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];		
	}
	std::vector<optype>::iterator op = base::opstr.begin();

	for(auto t=t_left.begin();t!=t_left.end();t++){
		if(op->o2 == -1){
			add_op_left(*op,*t);
		}
		else if(op->o1 > -2){
			remove_op_left(*op,*t);
		}
		else{
			base::sP[ op->o2 ] *= -1;
		}
		op++;
	}

	for(auto t=t_right.rbegin();t!=t_right.rend();t++){
		if(op->o2 == -1){
			add_op_right(*op,*t);
		}
		else if(op->o1 > -2){
			remove_op_right(*op,*t);
		}
		else{
			base::sP[ op->o2 ] *= -1;
		}
		op++;
	}
}


double im_base::update_times(int n,int n_update){
	n = std::min(n,std::min(M_left,M_right));
	n = std::max(1,n);
	int update=0;
	std::vector<double> times(n,0);
	for(int r=0;r<n_update;r++){
		int m = std::floor(base::ran()*(M_left-n+1));

		double t_1,t_2;
		if(m>0){t_1 = t_left[m-1];}
		else{t_1 = 0;}

		if(m+n<M_left){t_2 = t_left[m+n];}
		else{t_2 = t_f;}

		for(auto it=times.begin();it!=times.end();it++)
			*it=(t_1+base::ran()*(t_2-t_1));

		std::sort(times.begin(),times.end());

		double W_new=1;
		double W_old=1;
		
		std::vector<optype>::iterator op = base::opstr.begin()+m;
		std::vector<double>::iterator t1 = times.begin();
		std::vector<double>::iterator t2 = t_left.begin()+m;
		
		for(int p=0;p<n;p++){
			W_new *= weight(*op,*t1);
			W_old *= weight(*op,*t2);
			op++; t1++; t2++;
		}

		if(base::ran()*W_old<W_new){
			std::copy(times.begin(),times.end(),t_left.begin()+m);
			update++;
		}

		m = std::floor(base::ran()*(M_right-n+1));
		
		if(m>0){t_1 = t_right[m-1];}
		else{t_1 = 0;}

		if(m+n<M_right){t_2 = t_right[m+n];}
		else{t_2 = t_f;}
		
		for(auto it=times.begin();it!=times.end();it++)
			*it=(t_1+base::ran()*(t_2-t_1));

		std::sort(times.begin(),times.end());

		W_new=1;
		W_old=1;
		
		op = base::opstr.begin()+M_left+m;
		t1 = times.begin();
		t2 = t_right.begin()+m;

		for(int p=0;p<n;p++){
			W_new *= weight(*op,*t1);
			W_old *= weight(*op,*t2);
			op++; t1++; t2++;
		}

		if(base::ran()*W_old<W_new){
			std::copy(times.begin(),times.end(),t_right.begin()+m);
			update++;
		}
	}

	return update/double(2*n_update);
}




double inline im_base::time(int p){
	if(p<M_left){
		return t_left[p];
	}
	else{
		p -= M_left;
		return t_right[M_right-p-1];
	}
}



void im_base::midpoint(short_vec::iterator spins){
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
	int M_left_new = (Nop_left*4)/3;
	
	if(M_left_new > M_left){
		int nl = M_left_new - M_left;
		optype val;
		base::opstr.insert(base::opstr.begin(),nl,val);
		t_left.insert(t_left.begin(),nl,0);

		double t_min = t_left[0];
		std::vector<double>::iterator it=t_left.begin();
		for(int i=0;i<nl;i++){
			*it=base::ran()*t_min;
			it++;
		}

		std::sort(t_left.begin(),t_left.begin()+nl);
		M_left = M_left_new;
	}

	int M_right_new = (Nop_right*4)/3;
	if(M_right_new > M_right){	
		int nr = M_right_new - M_right;
		optype val;
		base::opstr.insert(base::opstr.end(),nr,val);
		t_right.insert(t_right.begin(),nr,0);

		double t_min = t_left[0];
		std::vector<double>::iterator it=t_right.begin();
		for(int i=0;i<nr;i++){
			*it=base::ran()*t_min;
			it++;
		}

		std::sort(t_right.begin(),t_right.begin()+nr);
		M_right = M_right_new;
	}

	base:M = base::opstr.size();
	X.resize(4*base::M);
}











#endif