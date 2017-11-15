#ifndef __opstr_iterator_INCLUDED__
#define __opstr_iterator_INCLUDED__

#include <algorithm>
#include <vector>
#include <tuple>

#include "optype.h"

typedef std::vector<signed char> Spin_vector;
typedef std::tuple<int,double,Spin_vector> Value;


// member typedefs provided through inheriting 0 std::iterator
class state_iterator: public std::iterator<std::input_iterator_tag,Value>{
	double msprop;
	int p;
	Spin_vector spins;
	std::vector<optype>::iterator opstr,opstr_end;

public:
	explicit state_iterator(Spin_vector spins_i,std::vector<optype>::iterator _opstr,std::vector<optype>::iterator _opstr_end) : opstr(_opstr), opstr_end(_opstr_end) {
		msprop = 0.0;
		p=0;
		for(auto s_i : spins_i){
			spins.push_back(s_i);
			msprop += s_i;
		}
		msprop /= 2.0;
	}
	state_iterator& operator++() {
		opstr++;
		p++;
		while(opstr->o1 > -2  && opstr != opstr_end){
			opstr++;
			// std::cout << opstr->o1 << "  " << opstr->o2 << std::endl;
			p++;
		}

		spins[opstr->o2] *= -1;
		msprop += spins[opstr->o2];

		return *this;
	}
	state_iterator operator++(int) {state_iterator retval = *this; ++(*this); return retval;}
	bool operator==(state_iterator other) const {return opstr == other.opstr;}
	bool operator!=(state_iterator other) const {return opstr != other.opstr;}
	Value operator*() const {
		return std::forward_as_tuple(p,msprop,spins);
	}
};




#endif