#ifndef __base_INCLUDED__
#define __base_INCLUDED__

#include <vector>
#include <stack>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
#include "state_iterator.h"
#include "optype.h"


class base{
	protected:
		int M;
		const int Fl;
		const int Fr;
		const int N;

		std::vector<signed char> sR;
		std::vector<signed char> sP;
		std::vector<signed char> sL;


		std::vector<int>  Vl;
		std::vector<int>  Vr;

		std::mt19937_64 gen;
		std::uniform_real_distribution<double> dist;

		std::vector<optype> opstr;
		std::stack<int> stk;
		std::vector<int> X;

		// private functions:
		void link_verticies();
		void visit_cluster();
		void flip_cluster();
		
		double inline ran(void);


	public:
		base(int,const int, const int,const int,
			const std::vector<signed char>,const std::vector<signed char>);
		base(int,const int);
		~base() {};
		std::vector<optype>::iterator opstr_begin() {return opstr.begin();}
		std::vector<optype>::iterator opstr_end() {return opstr.end();}
		state_iterator state_begin() {return state_iterator(sL,opstr.begin(),0,M);}
		state_iterator state_end() {return state_iterator(sL,opstr.end(),M+1,M);}
		void print_opstr(bool);
		void cluster_update();
		virtual void diagonal_update() = 0;
		int inline get_M(void) {return M;}
		int inline get_N(void)	{return N;}

};

base::base( int _M,
			const int _N
			) : M(_M), N(_N), Fl(0), Fr(0)
{

	//seeding random number generator
	unsigned int lo,hi,s;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	s=((unsigned long long)hi << 32) | lo;

	gen.seed(s);
	dist = std::uniform_real_distribution<double>(0.0,1.0);


	if((Fl==0) != (Fr==0)){
		std::cout << "Fr and Fl must both be equal to 0." << std::endl;
		exit(-1);
	}

	opstr.resize(M);
	X.resize(4*M);
	Vl.resize(N);
	Vr.resize(N);
	sP.resize(N);
	
	for(int i=0;i<N;i++){
		int s = 2*std::floor(2*ran())-1;
		sR.push_back(s);
	}

	for(auto s : sR){
		sL.push_back(s);
	}
}

base::base( int _M,
			const int _N,
			const int _Fl,
			const int _Fr,
			const std::vector<signed char> _sL,
			const std::vector<signed char> _sR
			) : M(_M), N(_N), Fl(_Fl), Fr(_Fr)
{

	//seeding random number generator
	unsigned int lo,hi,s;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	s=((unsigned long long)hi << 32) | lo;

	gen.seed(s);
	dist = std::uniform_real_distribution<double>(0.0,1.0);


	if((Fl==0) != (Fr==0)){
		std::cout << "Fr and Fl must both be equal to 0." << std::endl;
		exit(-1);
	}

	opstr.resize(M);
	X.resize(4*M);
	Vl.resize(N);
	Vr.resize(N);
	sP.resize(N);
	
	for(auto s : _sR){
		sR.push_back(s);
	}

	for(auto s : _sL){
		sL.push_back(s);
	}
}



double inline base::ran(void){
	return dist(gen);
}


void base::cluster_update(){
	link_verticies();

	// first mark clusters which are attached to end spins which can't be flipped
	// left side first

	if(Fl==-1){
		for(int i=0;i<N;i++){
			if(Vl[i]>=0){
				stk.push(Vl[i]);
				visit_cluster();
			}
		}
	}

	// right side next
	if(Fr==-1){
		for(int i=0;i<N;i++){
			if(Vr[i]>=0){
				stk.push(Vr[i]);
				visit_cluster();
			}
		}
	}

	// attempt with 50/50 probability to flip clusters which are attached to end spins which can be flipped
	// left side first
	if(Fl==1){
		for(int i=0;i<N;i++){
			if(Vl[i]>=0 && X[Vl[i]]>=0){
				stk.push(Vl[i]);
				if(ran()>=0.5){
					flip_cluster();
				}
				else{
					visit_cluster();
				}
			}
		}
	}

	// right side next
	if(Fr==1){
		for(int i=0;i<N;i++){
			if(Vr[i]>=0 && X[Vr[i]]>=0){
				stk.push(Vr[i]);
				if(ran()>=0.5){
					flip_cluster();
				}
				else{
					visit_cluster();
				}
			}
		}
	}

	// flipping bluk cluster with 50/50 probability on all clusters which have not been visited
	for(int v0=0;v0<4*M;v0++){
		if(X[v0]>=0){
			int p=v0/4;
			if(opstr[p].o1<0){
				int v1=X[v0];
				stk.push(v1);
				if(ran()>=0.5){
					opstr[p].o1^=1; // flip operator with some bit operation
					X[v0]=-2;
					flip_cluster();
				}
				else{
					X[v0]=-1;
					visit_cluster();
				}
			}
		}
	}

	// flipping spins at the end which a cluster hit.
	if(Fl==0 && Fr==0){
		for(int i=0;i<N;i++){
			int vl=Vl[i];
			if(vl==-1){
				if(ran()>=0.5){sL[i]*=-1;sR[i]*=-1;}
			}
			else{
				if(X[vl]==-2){sL[i]*=-1;sR[i]*=-1;}
			}
		}
	}
	else{
		for(int i=0;i<N;i++){
			int vl=Vl[i];	int vr=Vr[i];
			if(vl!=-1 && X[vl]==-2){sL[i]*=-1;}
			if(vr!=-1 && X[vr]==-2){sR[i]*=-1;}
			if(Fl==1 && Fr==1 && vr==-1 && vl==-1 && ran()>=0.5){ // flip spins disconnected from operators
				sL[i]*=-1;sR[i]*=-1;
			}

		}
	}
}



void base::visit_cluster(){
	while(!stk.empty()){
		int v=stk.top(); stk.pop();
		if(v >= 4*M) continue;
		X[v]=-1;
		int p=v/4;
		if(opstr[p].o1>=0){
			for(int v1=4*p; v1<4*p+4;v1++){
				if(X[v1]>=0){
					stk.push(X[v1]);
					X[v1]=-1;
				}// end if(X[v1]>=0)
			}// end for(int v1=v+1;...
		}//end if(opstr[p].o1>=0)
	}// end while(stk.empty()!) 
}



void base::flip_cluster(){
	while(!stk.empty()){
		int v=stk.top(); stk.pop();
		if(v >= 4*M) continue;
		X[v]=-2;
		int p=v/4;
		if(opstr[p].o1>=0){
			for(int v1=4*p; v1<4*p+4;v1++){
				if(X[v1]>=0){
					stk.push(X[v1]); 
					X[v1]=-2;
				}// end if(X[v1]>=0)
			}// end for(int v1=v+1;...
		}
		else{
			opstr[p].o1^=1; // flip operator with some bit operation
		}// end if(opstr[p].o1>=0)
	}// end while(stk.empty()!) 
}


// increasing p
// <Vl| -> |Vr>
void base::link_verticies(){
	for(int i=0; i<N; i++){ Vl[i] = Vr[i] = -1; }
	for(int i=0; i<4*M; i++){ X[i]=-1; }
	for(int p=0; p<M;p++){
		int o1=opstr[p].o1;
		int o2=opstr[p].o2;
		if(o2>=0){
			int v0=4*p;
			if(o1>-1){
				if(Vr[o1]==-1){ Vl[o1]=v0; }
				else{	X[Vr[o1]]=v0;
							X[v0]=Vr[o1];	}// end if(Vr[i0]==-1)
	
				if(Vr[o2]==-1){	Vl[o2]=v0+1; }
				else{	X[Vr[o2]]=v0+1;
						X[v0+1]=Vr[o2];	}// end if(Vr[i1]==-1)
				Vr[o2]=v0+3; 			Vr[o1]=v0+2;
			}
			else{
				if(Vr[o2]==-1){ Vl[o2]=v0; }
				else{ X[Vr[o2]]=v0;
							X[v0]=Vr[o2];	}// end if(Vr[o2]==-1)
				
				Vr[o2]=v0+2;
			}// end if(o1>=0)
		}//end if(o2>=0)
	}// end for(int p=0;...



	if(Fl!=0 && Fr!=0){
		for(int i=0;i<N;i++){
			if(Vl[i]>=0){ X[Vl[i]]=4*M+i;}
			if(Vr[i]>=0){ X[Vr[i]]=4*M+N+i;}
		}	
	}
	else if(Fl==0 && Fr==0){
		for(int i=0;i<N;i++){
			int f=Vl[i];
			if(f!=-1){ int l=Vr[i]; X[f]=l; X[l]=f;}
		}
	}
	else{ std::cout << "boundary mismatch" << std::endl; std::exit(11);}//end if(Fl!=0 && Fr!=0)
}


void base::print_opstr(bool link){
	std::cout << "p=   ";
	for(int p=0;p<M;p++){
		std::cout << p << " ";
	}
	std::cout << std::endl;
	for(int i=0;i<N;i++){
		std::cout << i << " ";
		std::cout << char(-sL[i]+44);
		std::cout << " ";
		for(int p=0;p<M;p++){
			std::cout << "-";
			int o1=opstr[p].o1;
			int o2=opstr[p].o2;
			if(o1 >= 0){
				if(o1==i || o2==i){ std::cout << "J";	}
				else{ std::cout << "-"; }
			}
			else{
				if(o1>-2){
					if(o2==i){ std::cout << "h";	}
					else{ std::cout << "-"; }
				}
				else{
					if(o2==i){ std::cout << "H";	}
					else{ std::cout << "-"; }
				}
			}
		}
		std::cout << "-";
		std::cout << " ";
		std::cout << char(-sR[i]+44) << " ";
		std::cout << i;
		std::cout << std::endl;
	}
	if(link){
		link_verticies();
		std::cout << std::setw(8) << std::endl;
		for(int p=M-1;p>=0;p--){
			for(int i=0;i<4;i++){
				std::cout << 4*p+i << std::setw(4) << X[4*p+i] << std::setw(8);
			}
			std::cout << std::endl;
		}
	}
}
#endif
