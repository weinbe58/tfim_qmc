#ifndef __QAQMC_INCLUDED__
#define __QAQMC_INCLUDED__

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


struct optype { public: int o1, o2; };


template<int L,int d>
struct get_N
{
	enum { N = get_N<L,d-1>::N * L};
};

template<int L>
struct get_N<L,1>
{
	enum { N = L };
};


template<int L,int d,int M,int mann>
class qaqmc{

	protected:
		enum{N=get_N<L,d>::N};
		enum{Nb=d*get_N<L,d>::N};

		double S_f,dS;
		const double dS_f;

		int Fl;
		int Fr;
		std::string output_file;

		double Eising,mi,ma,m2;

		signed char spinsR[N];
		signed char spins[N];
		signed char spinsL[N];

		int nns[Nb];
		int bst[2*Nb];
		int Vl[N];
		int Vr[N];

		//MTRand ran;

		std::mt19937_64 gen;
		std::uniform_real_distribution<double> dist;

		std::stack<int> stk;
		optype opstr[2*M];
		int X[8*M];

		// private functions:
		void diagonal_update();
		void inline update_S_f(int);
		void link_verticies();
		void cluster_update();
		void visit_cluster();
		void flip_cluster();
		void Measurement();

		void print_opstr(bool);

		double inline ran(void){
			return dist(gen);
		}


	public:
		qaqmc(const double,const std::string&);

		void write_out();
		void write_out_lock();
		void MCramp(int);

};


template<int L,int d,int M,int mann>
qaqmc<L,d,M,mann>::qaqmc(const double _S_f,
							  const std::string & _output_file
							 ) : dS_f(_S_f/mann)
{
	//seeding random number generator
	unsigned int lo,hi,s;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	s=((unsigned long long)hi << 32) | lo;

	gen.seed(s);
	dist = std::uniform_real_distribution<double>(0.0,1.0);

	output_file = _output_file;

	for(int i=0;i<N;i++){
		int LL = 1;
		for(int j=0;j<d;j++){
			int id=(i/LL)%L;
			int b = d*i+j;

			nns[b]=i+((id+1)%L-id)*LL;
			bst[2*b]=i;
			bst[2*b+1]=nns[b];
			LL*=L;
		}
	}

	for(int p=0;p<2*M;p++){ opstr[p].o1=-1; opstr[p].o2=int(std::floor(N*ran()));}

	/*	list of bc:
			px = polarized state in x direction
			pz = polarized state in z direction
			tr = trace 
	*/
	Fl=1; Fr=1;


}


template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::write_out_lock(){
	std::string filelock=output_file+".lock";
	std::ofstream lockstream;
	std::ofstream outstream; // Stream which the bin averages go to
	bool written=false;
	bool open=false;
	int count = 0;

	std::stringstream buffer;
	
	buffer << std::setprecision(14) << std::fixed;
	buffer << std::setw(20) << Eising;
	buffer << std::setw(20) << mi;
	buffer << std::setw(20) << ma;
	buffer << std::setw(20) << m2;
	buffer << std::endl;


 	while(!written){
		if(!std::ifstream(filelock.c_str())){ // if lock file exists wait continue loop.
			lockstream.open(filelock.c_str()); // open lock file to indicate that file is open and no other program can write.
			while(!open){ 
				outstream.open(output_file.c_str(), std::ios::app);
				open=outstream.is_open();
				if(open){	
					outstream << buffer.str() << std::flush;
					outstream.close();
					written=true;
					lockstream.close();
					int del=std::remove(filelock.c_str());
					while(del != 0){del=std::remove(filelock.c_str());} 					
				}
				count++;
				if(count>10000){std::cout << "maximum number of attemps to write to file reached!" << std::endl; exit(4);}
			}
		}
		count++;
		if(count>10000){std::cout << "maximum number of attemps to write to file reached!" << std::endl; exit(4);}
	}
}


template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::write_out(){
	std::ofstream outstream; // Stream which the bin averages go to
	bool written=false;
	bool open=false;
	int count=0;
	std::stringstream buffer;
	
	buffer << std::setprecision(14) << std::fixed;
	buffer << std::setw(20) << Eising;
	buffer << std::setw(20) << mi;
	buffer << std::setw(20) << ma;
	buffer << std::setw(20) << m2;
	buffer << std::endl;


	while(!open){ 
		outstream.open(output_file.c_str(), std::ios::app);
		open=outstream.is_open();
		if(open){	
			outstream << buffer.str() << std::flush;
			outstream.close();
			written=true;					
		}
		count++;
		if(count>10000){std::cout << "maximum number of attemps to write to file reached!" << std::endl; exit(4);}
	}
}


template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::MCramp(const int npbin){
	
	for(int k=0;k<npbin;k++){
		std::cout << k << std::endl;
		Eising=0;
		mi=0;
		ma=0;
		m2=0;
		for(int p=0;p<2*M;p++){ opstr[p].o1=-1; opstr[p].o2=int(std::floor(N*ran()));}
		for(int i=0;i<N;i++){
			spinsL[i]=2*int(std::floor(ran()))-1;
			spinsR[i]=2*int(std::floor(ran()))-1;
		}
		S_f = 0.0;
		dS = 0.0;
		for(int j=0;j<mann;j++){
			diagonal_update();
			cluster_update();
			update_S_f(j);
		}
		Measurement();

		Eising /= (N*npbin);
		mi /= npbin;
		ma /= npbin;
		m2 /= npbin;
	}
}




template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::diagonal_update(){
	for(int i=0;i<N;i++){
		spins[i]=spinsL[i];
	}// end for(int i=1;...
	double S=0;
	for(int p=0;p<M;p++){
		if(opstr[p].o1>-2){
			bool accept=false;
			S += dS;
			double p1 = (1-S);
			double p2 = ((2*d-1)*S+1);
			while(true){
				if(ran()*p2<p1){
					opstr[p].o1=-1; 
					opstr[p].o2=int(ran()*N);
					break;
				}
				else{
					int ib = 2*int(ran()*Nb);
					int i = bst[ib];
					int j = bst[ib+1];
					if(spins[i]*spins[j]>0){
						opstr[p].o1=i;
						opstr[p].o2=j;
					}
					break;
				}
			}
		}
		else{
			spins[ opstr[p].o2 ] *= -1;
		}
	}
	for(int p=M;p<2*M;p++){
		if(opstr[p].o1>-2){
			bool accept=false;
			double p1 = (1-S);
			double p2 = ((2*d-1)*S+1);
			while(true){
				if(ran()*p2<p1){
					opstr[p].o1=-1; 
					opstr[p].o2=int(ran()*N);
					break;
				}
				else{
					int ib = 2*int(ran()*Nb);
					int i = bst[ib];
					int j = bst[ib+1];
					if(spins[i]*spins[j]>0){
						opstr[p].o1=i;
						opstr[p].o2=j;
					}
					break;
				}
			}
		}
		else{
			spins[ opstr[p].o2 ] *= -1;
		}
		S -= dS;
	}
}




template<int L,int d,int M,int mann>
void inline qaqmc<L,d,M,mann>::update_S_f(int k){
	S_f += dS_f;
	dS += dS_f/M;
}


template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::Measurement(){
	double mprop;
	for(int i=0;i<N;i++){spins[i]=spinsL[i];}
	mprop=0;
	for(int i=0;i<N;i++){mprop+=spins[i];}// end for(int i=0;i<N;i++)

	mprop=mprop/2;
	
	for(int p=0;p<M;p++){
		if(opstr[p].o1<-1){
			spins[opstr[p].o2]*=-1;
			mprop+=spins[opstr[p].o2];
		}
	}

	double ising=0;
	for(int b=0;b<Nb;b++){
		int i=bst[2*b];
		int j=bst[2*b+1];
		ising += spins[i]*spins[j];
	}
	Eising += ising;

	double Mag=mprop/N;
	mi+=Mag;
	ma+=std::abs(Mag);
	m2+=Mag*Mag;

}


template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::cluster_update(){
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

	// flipping cluster with 50/50 probability on all clusters which have not been visited
	for(int v0=0;v0<8*M;v0++){
		if(X[v0]>=0){
			int p=v0/4;
			if(opstr[p].o1<0){
				int v1=X[v0];
				stk.push(v1);;
				if(ran()>=0.5){
					opstr[p].o1=(opstr[p].o1^1); // flip operator with some bit operation
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
				if(ran()<0.5){spinsL[i]*=-1;spinsR[i]*=-1;}
			}
			else{
				if(X[vl]==-2){spinsL[i]*=-1;spinsR[i]*=-1;}
			}
		}
	}
	else{
		for(int i=0;i<N;i++){
			int vl=Vl[i];	int vr=Vr[i];
			if(vl!=-1){ if(X[vl]==-2){spinsL[i]*=-1;} }
			if(vr!=-1){ if(X[vr]==-2){spinsR[i]*=-1;} }
		}
	}
}


template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::visit_cluster(){
	while(!stk.empty()){
		int v=stk.top(); stk.pop();
		if(v >= 8*M) continue;
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


template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::flip_cluster(){
	while(!stk.empty()){
		int v=stk.top(); stk.pop();
		if(v >= 8*M) continue;
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
			opstr[p].o1=(opstr[p].o1^1); // flip operator with some bit operation
		}// end if(opstr[p].o1>=0)
	}// end while(stk.empty()!) 
}


// increasing p
// <Vl| -> |Vr>
template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::link_verticies(){
	for(int i=0; i<N; i++){ Vl[i] = Vr[i] = -1; }
	for(int i=0; i<8*M; i++){ X[i]=-1; }
	for(int p=0; p<2*M;p++){
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
			if(Vl[i]>=0){ X[Vl[i]]=8*M+i;}
			if(Vr[i]>=0){ X[Vr[i]]=8*M+N+i;}
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


template<int L,int d,int M,int mann>
void qaqmc<L,d,M,mann>::print_opstr(bool link){
	std::cout << "p=   ";
	for(int p=0;p<2*M;p++){
		std::cout << p << " ";
	}
	std::cout << std::endl;
	for(int i=0;i<N;i++){
		std::cout << i << " ";
		std::cout << char(-spinsL[i]+44);
		std::cout << " ";
		for(int p=0;p<2*M;p++){
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
		std::cout << char(-spinsR[i]+44) << " ";
		std::cout << i;
		std::cout << std::endl;
	}
	if(link){
		link_verticies();
		std::cout << std::setw(8) << std::endl;
		for(int p=2*M-1;p>=0;p--){
			for(int i=0;i<4;i++){
				std::cout << 4*p+i << std::setw(4) << X[4*p+i] << std::setw(8);
			}
			std::cout << std::endl;
		}
	}
}

#endif
