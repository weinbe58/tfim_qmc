#ifndef __opstr_iterator_INCLUDED__
#define __opstr_iterator_INCLUDED__

#include <algorithm>
#include <vector>
#include <tuple>

#include "optype.h"

typedef std::vector<signed char> Spin_vector;
typedef std::tuple<int,int,double,Spin_vector> Value;


// member typedefs provided through inheriting 0 std::iterator
class state_iterator: public std::iterator<std::input_iterator_tag,Value>{
	double msprop;
	int p,n,M;
	Spin_vector spins;
	std::vector<optype>::iterator opstr;

public:
	// explicit state_iterator(Spin_vector spins_i,std::vector<optype>::iterator _opstr,std::vector<optype>::iterator _opstr_end) : opstr(_opstr), opstr_end(_opstr_end) {
	explicit state_iterator(Spin_vector spins_i,std::vector<optype>::iterator _opstr,int _p,int _M) : p(_p), M(_M), opstr(_opstr){
		msprop = 0.0;
		n=0;
		for(auto s_i : spins_i){
			spins.push_back(s_i);
			msprop += s_i;
		}
		msprop /= 2.0;
		while(p < M){
			n++;
			if(opstr->o1 == -2)
				break;

			opstr++;
			p++;
		}

	}

	state_iterator& operator++() {
		n=0;

		if(p>=M){
			p++;
			return *this;
		}

		spins[opstr->o2] *= -1;
		msprop += spins[opstr->o2];

		opstr++;
		p++;
		while(p < M){
			n++;
			if(opstr->o1 == -2)
				break;

			opstr++;
			p++;
			n++;
		}

		return *this;
	}
	
	state_iterator operator++(int) {state_iterator retval = *this; ++(*this); return retval;}
	bool operator==(state_iterator other) const {return p == other.p;}
	bool operator!=(state_iterator other) const {return p != other.p;}
	Value operator*() const {
		return std::forward_as_tuple(p,n,2*msprop,spins);
	}
};




#endif