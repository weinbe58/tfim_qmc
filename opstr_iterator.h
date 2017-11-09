#ifndef __opstr_iterator_INCLUDED__
#define __opstr_iterator_INCLUDED__

#include <algorithm>
#include <vector>
#include <tuple>

struct optype { public: int o1=0; int o2=0; };

typedef std::vector<signed char> Spin_vector;
typedef std::tuple<optype,double,Spin_vector> Value;


// member typedefs provided through inheriting 0 std::iterator
class opstr_iterator: public std::iterator<std::input_iterator_tag,Value>{
	double msprop;
	Spin_vector spins;
	std::vector<optype>::iterator opstr;

public:
	explicit opstr_iterator(Spin_vector spins_i,std::vector<optype>::iterator _opstr) : opstr(_opstr) {
		msprop = 0.0;
		for(auto s_i : spins_i){
			spins.push_back(s_i);
			msprop += s_i;
		}
		msprop /= 2.0;
	}
	iterator& operator++() {
		++opstr;
		int o1 = opstr->o1;
		int o2 = opstr->o2;
		if(o1 < -1){
			spins[o2] *= -1;
			msprop += spins[o2];
		}

		return *this;
	}
	iterator operator++(int) {iterator retval = *this; ++(*this); return retval;}
	bool operator==(iterator other) const {return opstr == other.opstr;}
	bool operator!=(iterator other) const {return !(*this == other);}
	Value operator*() const {
		optype op = *opstr;		
		return std::forward_as_tuple(op,msprop,spins);
	}
};

#endif