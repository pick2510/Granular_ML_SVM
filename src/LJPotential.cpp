#include "LJPotential.h"
#include <cmath>
#include <cstring>
	LJPotential::~LJPotential()
	{
	}
	double LJPotential::IsotropicEnergyFormula(double distance, const AtomInfo & a1, const AtomInfo & a2)
	{
		double a=this->sigma*this->sigma/distance/distance;
		a=a*a*a;
		return this->epsilon4*(a*a-a)-this->Shift;
	}
	LJPotential::LJPotential(DimensionType Dimension, double epsilon, double sigma, double Rcut):IsotropicPairPotential(Dimension, Rcut)
	{
		this->epsilon4=4*epsilon;
		this->sigma=sigma;
		this->Shift=0.0;
		this->Shift=this->IsotropicEnergyFormula(Rcut, "", "");
	}
	double LJPotential::IsotropicForceFormula(double distance, const AtomInfo & a1, const AtomInfo & a2)
	{
		double temp = this->sigma / distance;
		return (-1.0)*this->epsilon4 * (6 * std::pow(temp, (double)(7)) - 12 * std::pow(temp, (double)(13))) / this->sigma;
	}
	double LJPotential::IsotropicSecondDerivative(double distance, const AtomInfo & a1, const AtomInfo & a2)
	{
		double temp = this->sigma / distance;
		return (-1.0)*this->epsilon4 * (42 * std::pow(temp, (double)(8)) - 12 * 13 * std::pow(temp, (double)(14))) / this->sigma / this->sigma;
	}
