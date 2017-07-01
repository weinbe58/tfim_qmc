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


template<int L,int d,int M,int Nm,int mstep>
class qaqmc{

	enum{N=get_N<L,d>::N};
	enum{Nb=d*get_N<L,d>::N};

	const double S_i;
	const double S_f;
	const double dS;
	int Fl;
	int Fr;
	double S;
	std::string output_file;

	int Ef[Nm+1];
	int Et[Nm+1];
	double Eising[Nm+1];
	double E[Nm+1];
	double mi[Nm+1];
	double ma[Nm+1];
	double m2[Nm+1];

	signed char spinsR[N];
	signed char spins[N];
	signed char spinsL[N];

	int nns[Nb];
	int bst[2*Nb];
	int p_list[Nm+1];
	int Vl[N];
	int Vr[N];

	//MTRand ran;

	std::mt19937_64 gen;
	std::uniform_real_distribution<double> dist;


	std::vector<optype> opstr;
	std::stack<int> stk;
	std::vector<int> X;
	std::vector<bool> Measure;

	// private functions:
	void diagonal_update();
	void inline update_S(int);
	void link_verticies();
	void cluster_update();
	void visit_cluster();
	void flip_cluster();

	void beginMeasure();
	void Measurement();
	void endMeasure();
	void double_opstr();
	void print_opstr(bool);

	double inline ran(void){
		return dist(gen);
	}


	public:
		qaqmc(const double,const double,const std::string&,const std::string&);

		void write_out();
		void write_out_lock();
		void MCstep();
		void EQstep();

};


template<int L,int d,int M,int Nm,int mstep>
qaqmc<L,d,M,Nm,mstep>::qaqmc(const double _S_i,
							  const double _S_f,
							  const std::string & bc,
							  const std::string & _output_file
							 ) : S_f(_S_f), S_i(_S_i), dS((_S_f-_S_i)/M)
{
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

	opstr.resize(2*M);
	Measure.resize(2*M+1);
	X.resize(8*M);

	int dm=0;

	if(Nm>2*M){ std::cout << "Nm is too large." << std::endl; std::exit(3);}
	else if(Nm<=0){ dm=0; }
	else{ dm=2*M/Nm; }


	for(int i=0;i<=2*M;i++){ Measure[i]=false; }

	int k=0;
	for(int i= -Nm/2;i<=Nm/2;i++){ 
		Measure[M+i*dm]=true;
		p_list[k]=M+i*dm;
		k++;
	}

	for(int p=0;p<2*M;p++){ opstr[p].o1=-1; opstr[p].o2=int(std::floor(N*ran()));}

	/*	list of bc:
			px = polarized state in x direction
			pz = polarized state in z direction
			tr = trace 
	*/
	if(bc=="pzpz"){
		Fl=-1; Fr=-1;
		for(int i=0;i<N;i++){
			spinsL[i]=1;
			spinsR[i]=1;
		}
	}
	else if(bc=="pxpx"){
		Fl=1; Fr=1;
		for(int i=0;i<N;i++){
			spinsL[i]=1;
			spinsR[i]=1;
		}
	}
	else if(bc=="tr"){
		Fl=0; Fr=0;
		for(int i=0;i<N;i++){
			spinsL[i]=1;
			spinsR[i]=1;
		}
	}
	else{
		std::cout << "init_states error: bc not recognized." << std::endl;
		std::exit(2);
	}// end bc if tree

	//seeding random number generator
	unsigned int lo,hi,s;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	s=((unsigned long long)hi << 32) | lo;

	gen.seed(s);
	dist = std::uniform_real_distribution<double>(0.0,1.0);
}


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::write_out_lock(){
	std::string filelock=output_file+".lock";
	std::ofstream lockstream;
	std::ofstream outstream; // Stream which the bin averages go to
	bool written=false;
	bool open=false;
	int count = 0;

	std::stringstream buffer;
	buffer << std::setprecision(14) << std::fixed;

	for(int i=0;i<Nm+1;i++){ 
		update_S(p_list[i]);
		buffer << std::setw(20) << S;
		buffer << std::setw(20) << E[i];
		buffer << std::setw(20) << mi[i];
		buffer << std::setw(20) << ma[i];
		buffer << std::setw(20) << m2[i];
	}
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


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::write_out(){
	std::ofstream outstream; // Stream which the bin averages go to
	bool written=false;
	bool open=false;
	int count=0;
	std::stringstream buffer;
	
	buffer << std::setprecision(14) << std::fixed;
	for(int i=0;i<Nm+1;i++){ 
		update_S(p_list[i]);
		buffer << std::setw(20) << S;
		buffer << std::setw(20) << E[i];
		buffer << std::setw(20) << mi[i];
		buffer << std::setw(20) << ma[i];
		buffer << std::setw(20) << m2[i];
	}
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


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::EQstep(){
	for(int j=0;j<mstep;j++){
		diagonal_update();
		cluster_update();
	}
}


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::MCstep(){
	beginMeasure();
	for(int j=0;j<mstep;j++){
		diagonal_update();
		cluster_update();
		Measurement();
	}
	endMeasure();
}


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::diagonal_update(){
	for(int i=0;i<N;i++){
		spins[i]=spinsL[i];
	}// end for(int i=1;...

	for(int p=0;p<2*M;p++){
		if(opstr[p].o1>-2){
			bool accept=false;
			update_S(p);
			double p1 = (1-S);
			double p2 = ((2*d-1)*S+1);
			while(!accept){
				if(ran()*p2<p1){
					opstr[p].o1=-1; 
					opstr[p].o2=int(ran()*N);
					accept=true;
				}
				else{
					int ib = 2*int(ran()*Nb);
					int i = bst[ib];
					int j = bst[ib+1];
					if(spins[i]*spins[j]>0){
						opstr[p].o1=i;
						opstr[p].o2=j;
					}
					accept=true;
				}
			}
		}
		else{
			spins[ opstr[p].o2 ] *= -1;
		}
	}
}


template<int L,int d,int M,int Nm,int mstep>
void inline qaqmc<L,d,M,Nm,mstep>::update_S(int p){
	S=(p<M ? S_i+dS*(p+1) : S_f-dS*(p-M));
}


template<int L,int d,int M, int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::double_opstr(){
	std::vector<optype> opstr_1;
	opstr_1.resize(2*M);
	for(int p=0;p<2*M;p++){opstr_1[p]=opstr[p];}
	opstr.resize(4*M);
	for(int p=0;p<2*M;p++){
		if(opstr_1[p].o1>-2){
			opstr[2*p]=opstr_1[p];
			opstr[2*p+1]=opstr_1[p];
		}
		else{
			opstr[2*p].o1=-1;
			opstr[2*p].o2=opstr_1[p].o2;
			opstr[2*p+1]=opstr_1[p];
		}
	}
}


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::Measurement(){
	double mprop;
	int k;
	for(int i=0;i<N;i++){spins[i]=spinsL[i];}
	mprop=0;
	for(int i=0;i<N;i++){mprop+=spins[i];}// end for(int i=0;i<N;i++)

	mprop=mprop/2;
	k=0;
	
	int o1,o2;
	for(int p=0;p<2*M;p++){
		o1=opstr[p].o1;
		o2=opstr[p].o2;
		if(Measure[p]){
			int ising = 0;
			for(int b=0;b<Nb;b++){
				int i=bst[2*b];
				int j=bst[2*b+1];
				ising += spins[i]*spins[j];
			}
			Eising[k] += ising;

			if(o1==-1){Et[k]++;}
			else if(o1==-2){Ef[k]++;}

			double Mag=mprop/N;
			mi[k]+=Mag;
			ma[k]+=std::abs(Mag);
			m2[k]+=Mag*Mag;
			k++;
		}
		if(o1<-1){
			spins[o2]*=-1;
			mprop+=spins[o2];
		}// end if(opstr[p].o1==-2)
	}// for(int p=0;p<2*M;p++)


	if(Measure[2*M]){
		int ising = 0;
		for(int b=0;b<Nb;b++){
			int i=bst[2*b];
			int j=bst[2*b+1];
			ising+=spins[i]*spins[j];
		}
		Eising[k] += ising;

		if(o1==-1){Et[k]++;}
		else if(o1==-2){Ef[k]++;}

		double Mag=mprop/N;
		mi[k]+=Mag;
		ma[k]+=std::abs(Mag);
		m2[k]+=Mag*Mag;
	}
}


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::beginMeasure(){
	for(int i=0;i<=Nm;i++){
		Eising[i]=0.0;
		Ef[i]=0.0;
		Et[i]=0.0;
		mi[i]=0.0;
		ma[i]=0.0;
		m2[i]=0.0;
	}
}


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::endMeasure(){
	for(int i=0;i<=Nm;i++){
		update_S(p_list[i]);
		
		Eising[i]/=double(N*mstep);
		E[i]= -S*Eising[i] - (1.0-S)*double(Ef[i])/(double(Et[i])+1.1e-16);

		mi[i]=mi[i]/mstep;
		ma[i]=ma[i]/mstep;
		m2[i]=m2[i]/mstep;
	}
}


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::cluster_update(){
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


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::visit_cluster(){
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


template<int L,int d,int M, int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::flip_cluster(){
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
template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::link_verticies(){
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


template<int L,int d,int M,int Nm,int mstep>
void qaqmc<L,d,M,Nm,mstep>::print_opstr(bool link){
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
