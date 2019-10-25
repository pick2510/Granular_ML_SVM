#ifndef LJPOTENTIAL_INCLUDED
#define LJPOTENTIAL_INCLUDED

#include "Potential.h"

class LJPotential: public IsotropicPairPotential
{
public:
	double epsilon4;//4*epsilon
	double sigma;
	double Shift;
	virtual ~LJPotential();
	LJPotential(DimensionType Dimension, double epsilon, double sigma, double Rcut);
	virtual double IsotropicEnergyFormula(double distance, const AtomInfo & a1, const AtomInfo & a2);
	virtual double IsotropicForceFormula(double distance, const AtomInfo & a1, const AtomInfo & a2);
	virtual double IsotropicSecondDerivative(double distance, const AtomInfo & a1, const AtomInfo & a2);
};

#endif