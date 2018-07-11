#ifndef __uniform_dist_INCLUDED__
#define __uniform_dist_INCLUDED__

#include <random>


class uniform_dist
{
	std::mt19937_64 gen;
	std::uniform_real_distribution<double> dist;

	public:


		uniform_dist(){
			//seeding random number generator
			unsigned int lo,hi,s;
			__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
			s=((unsigned long long)hi << 32) | lo;

			gen.seed(s);
			dist = std::uniform_real_distribution<double>(0.0,1.0);
		};

		uniform_dist(double low, double high){
			//seeding random number generator
			unsigned int lo,hi,s;
			__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
			s=((unsigned long long)hi << 32) | lo;

			gen.seed(s);
			dist = std::uniform_real_distribution<double>(low,high);
		};

		uniform_dist(double low, double high,unsigned int s){
			gen.seed(s);
			dist = std::uniform_real_distribution<double>(low,high);
		};

		inline double operator()(void){
			return dist(gen);
		};
};



#endif