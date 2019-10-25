#ifndef RANDOMGENERATOR_INCLUDED
#define RANDOMGENERATOR_INCLUDED

#include <random>

class RandomGenerator
{
private:
	std::mt19937 gen;
	std::uniform_real_distribution<double> dist;
	void initialize(void);
public:
	void seed(int s);
	RandomGenerator();
	RandomGenerator(int Seed);
	double RandomDouble(void);
};
void GetRandomVector(int Dimension, RandomGenerator & gen, double * result);

#endif
