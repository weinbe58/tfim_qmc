
#include "qaqmc.h"

template<int L,int d,int M,int Nm,int mstep>
class qmc : qaqmc<L,d,M,Nm,mstep>{

	protected:
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


	public:
		qmc(const double,const std::string&,const std::string&);

};

template<int L,int d,int M,int Nm,int mstep>
qmc<L,d,M,Nm,mstep>::qmc(const double _S,
						 const std::string & bc,
						 const std::string & _output_file
						 ) : S_f(_S), S_i(_S), dS(0)
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

	for(int i=0;i<=2*M;i++){ Measure[i]=false; }

	int k=0;
	for(int i= -Nm/2;i<=Nm/2;i++){ 
		Measure[M+i]=true;
		p_list[k]=M+i;
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


