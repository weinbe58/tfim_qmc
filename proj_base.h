#ifndef __proj_base_INCLUDED__
#define __proj_base_INCLUDED__

#include "base.h"
#include <algorithm>
#include <vector>

class proj_base : public base{
	protected:
		virtual void move_op(int) = 0;

	public:
		proj_base(int,const int,const int, const int, const std::vector<signed char>,const std::vector<signed char>);
		proj_base(int,const int, const std::vector<signed char>,const std::vector<signed char>,
				  const std::vector<signed char>,const std::vector<signed char>);
		proj_base(int,const int,const int, const int);
		proj_base(int,const int);
		~proj_base() {};

		void diagonal_update(void);
		void diagonal_update(std::vector<signed char>::iterator);
		void diagonal_update(std::vector<signed char>::iterator,int&,int&,int&);
		void double_M();
		void midpoint(std::vector<signed char>::iterator);
		void initialize_opstr();
};


proj_base::proj_base( int _M,
			const int _N,
			const std::vector<signed char> _Fl,
			const std::vector<signed char> _Fl,
			const std::vector<signed char> _sL,
			const std::vector<signed char> _sR) : base::base(_M,_N,_Fl,_Fr,_sL,_sR) {
	initialize_opstr();
}

proj_base::proj_base( int _M,
			const int _N,
			const int _Fl,
			const int _Fr,
			const std::vector<signed char> _sL,
			const std::vector<signed char> _sR) : base::base(_M,_N,_Fl,_Fr,_sL,_sR) {
	initialize_opstr();
}

proj_base::proj_base( int _M,
			const int _N,
			const int _Fl,
			const int _Fr) : base::base(_M,_N,_Fl,_Fr) {
	initialize_opstr();
}

proj_base::proj_base( int _M, const int _N) : base::base(_M,_N) {
	initialize_opstr();
}


void proj_base::initialize_opstr(){
	for(int p=0;p<base::M;p++){
		base::opstr[p].o1=-1;
		base::opstr[p].o2=std::floor(base::N*base::ran());
	}
}

void proj_base::midpoint(std::vector<signed char>::iterator spins){
	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];
	}

	int p = 0;
	auto op = base::opstr.begin();

	while(p<=base::M/2){
		if(op->o1==-2){
			base::sP[ op->o2 ] *= -1;
		}
		op++;
		p++;
	}

	std::copy(base::sP.begin(),base::sP.end(),spins);
	
}

void proj_base::diagonal_update(){
	// this->print_opstr(false);
	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];
		// std::cout << int(base::sL[i]) << "  ";
	}// end for(int i=1;...
	// std::cout << std::endl;

	for(int p=0;p<base::M;p++){
		if(base::opstr[p].o1>-2){
			move_op(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}// end (opstr_l[p].o2<0)
	}
}

void proj_base::diagonal_update(std::vector<signed char>::iterator spins){
	// this->print_opstr(false);
	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];
		// std::cout << int(base::sL[i]) << "  ";
	}// end for(int i=1;...
	// std::cout << std::endl;
	int MM = base::M/2;
	for(int p=0;p<MM;p++){
		if(base::opstr[p].o1>-2){
			move_op(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}// end (opstr_l[p].o2<0)
	}
	std::copy(sP.begin(),sP.end(),spins);
	for(int p=MM;p<base::M;p++){
		if(base::opstr[p].o1>-2){
			move_op(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}// end (opstr_l[p].o2<0)
	}
}


void proj_base::diagonal_update(std::vector<signed char>::iterator spins,int &np,int &nm,int &nt){
	// this->print_opstr(false);
	for(int i=0;i<base::N;i++){
		base::sP[i]=base::sL[i];
		// std::cout << int(base::sL[i]) << "  ";
	}// end for(int i=1;...
	// std::cout << std::endl;
	int MM = base::M/2;
	for(int p=0;p<MM;p++){
		if(base::opstr[p].o1>-2){
			move_op(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}// end (opstr_l[p].o2<0)
	}
	std::copy(sP.begin(),sP.end(),spins);
	if(base::opstr[MM].o1==-2){
		if(base::sP[ base::opstr[MM].o2 ]==-1){np++;}
		else{nm++;}
	}
	else if(base::opstr[MM].o1==-1){
		nt++;
	}
	for(int p=MM;p<base::M;p++){
		if(base::opstr[p].o1>-2){
			move_op(p);
		}
		else{
			base::sP[ base::opstr[p].o2 ] *= -1;
		}// end (opstr_l[p].o2<0)
	}
}

void proj_base::double_M(){
	int Mnew = 2*base::M;
	std::vector<optype> opstr_temp(base::opstr);
	base::opstr.resize(Mnew);
	for(int p=0;p<base::M;p++){
		int o1 = opstr_temp[p].o1;
		int o2 = opstr_temp[p].o2;
		if(o1>-1){
			base::opstr[2*p].o1 = opstr_temp[p].o1;
			base::opstr[2*p].o2 = opstr_temp[p].o2;
			base::opstr[2*p+1].o1 = opstr_temp[p].o1;
			base::opstr[2*p+1].o2 = opstr_temp[p].o2;
		}
		else{
			base::opstr[2*p].o1 = opstr_temp[p].o1;
			base::opstr[2*p].o2 = opstr_temp[p].o2;
			base::opstr[2*p+1].o1 = -1;
			base::opstr[2*p+1].o2 = opstr_temp[p].o2;			
		}
	}
	base::M = Mnew;
}

#endif
