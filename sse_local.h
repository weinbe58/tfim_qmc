#ifndef __sse_local_INCLUDED__
#define __sse_local_INCLUDED__


class sse_local : sse_base
{
	protected:
		const int * bst;
		const int Nb;
		const double J,h,WJ,Wh;
		double W,beta;

		bool add_op(int);
		bool remove_op(int);

	public:
		sse_local(const int, const int, const * int,const double, const double, const double);
		~sse_local();

};

sse_local::sse_local(const int _N,const int _Nb, const int * bst, const double _J,const double _h, const double _beta) 
: sse_base::sse_base(_N), Nb(_Nb), bst(_bst), J(_J), h(_h), beta(_beta), WJ(_Nb*std::abs(J)), Wh(_N*abs(h)) {
	W = beta * ( WJ + Wh);
}


sse_local::add_op(int p){
	if(base::ran()*(base::M - sse_base::Nop + W) < W){
		while(true){
			if(ran()*WJ<Wh){
				base::opstr[p].o1=-1; 
				base::opstr[p].o2=int(base::ran()*base::N);
				break;
			}
			else{
				int ib = 2*int(base::ran()*Nb);
				int i = bst[ib];
				int j = bst[ib+1];
				if(spins[i]*spins[j]>0){
					base::opstr[p].o1=i;
					base::opstr[p].o2=j;
					break;
				}
			}
		}
		sse_base::Nop++;
	}
}

sse_local::remove_op(){
	if(ran()*(M - Nop + 1 + W) < (M - Nop + 1)){
		Nop--;
		base::opstr[p].o1=-1;
		base::opstr[p].o2=-1;
	}
}



#endif