#ifndef STRUCTUREOPTIMIZATION_INCLUDED
#define STRUCTUREOPTIMIZATION_INCLUDED

#include "PeriodicCellList.h"
#include "Potential.h"
#include <cassert>

//default: largest size_t
//if not zero, then partially optimized configurations are saved every these steps to StructureOptimization_Progress.ConfigPack
extern size_t StructureOptimization_SaveConfigurationInterval;
//default: false
extern bool StructureOptimization_SaveNumberedConfiguration;

extern bool StructureOptimization_SaveEquiEnergyConfiguration;
extern double StructureOptimization_SaveEquiEnergyConfigurationEnergyInterval;
extern std::string StructureOptimization_SaveEquiEnergyConfigurationName;


class NotANumberFound : public std::exception
{
public:
	NotANumberFound()
	{
	}
	virtual const char * what() const throw()
	{
		return "Caught a NaN!\n";
	}
};

//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
void RelaxStructure_MINOP_withoutPertubation(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=20000, double Emin=-1e10);

//do MINOP relax without pertubation, then with pertubation, return the best one
void RelaxStructure_MINOP(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=20000);

//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
void RelaxStructure_NLOPT(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);
//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
void RelaxStructure_NLOPT_Emin(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double Emin);
//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
void RelaxStructure_SteepestDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000);
void RelaxStructure_ConjugateGradient(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep=10000, double minG=1e-13, double minF=(-1.0)*::MaxEnergy);

//this function returns the physical "time," i.e., the sum of (step size divided by the modulus of the gradient)
double RelaxStructure_LocalGradientDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep = 10000, double maxError = 0.01, double defaultStepSize=1e-3, double maxTime = std::numeric_limits<double>::max() );

void RelaxStructure_LocalGradientDescent_StepSnapshots(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double maxError, double stepInterval, std::string packName);//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both


void RelaxStructure_FIRE(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep = 10000, double dtMax = 10);
void RelaxStructure_FIRE_Emin(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double dtMax, double Emin);//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both

//slower but does not need derivative information
void RelaxStructure_NLOPT_NoDerivative(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, double rescale=1.0, size_t MaxStep=100000, double fmin=(-1.0)*MaxEnergy);

inline void RelaxStructure(Configuration & List, Potential & pot, double Pressure, double MinDistance)
{
	try
	{
		RelaxStructure_NLOPT(List, pot, Pressure, 2, MinDistance);
	}
	catch(EnergyDerivativeToBasisVectors_NotImplemented a)
	{
		for(int i=0; i<5; i++)
		{
			RelaxStructure_NLOPT(List, pot, Pressure, 0, MinDistance);
			RelaxStructure_NLOPT_NoDerivative(List, pot, Pressure, 1, MinDistance);
		}
		RelaxStructure_NLOPT_NoDerivative(List, pot, Pressure, 2, MinDistance);
	}
}

//calculate Elastic Constants
//relax the structure before calling this function
//\epsilon_{ij}=C_{ijkl}*e_{kl}
double ElasticConstant(Potential & pot, Configuration structure, DimensionType i, DimensionType j, DimensionType k, DimensionType l, double Pressure, bool InfiniteTime=false);

//if presult is not nullptr, then presult[0..dimension^4] will be filled with C_{ijkl} 
void PrintAllElasticConstants(std::ostream & out, Configuration stru, Potential & pot, double Pressure, bool InfiniteTime = false, double * presult=nullptr);


//optimize the parameters for a elastic model and rotations so that it fits the target
void ElasticOptimize(std::function<BigVector(const std::vector<double> & param)> getElasticFunc, int nparam, std::vector<double> & param, std::vector<double> ubound, std::vector<double> lbound, const Configuration & c, Potential & pot, double pressure);


//if pEigenVectors is not nullptr, then fill it with eigenvectors
void GetHessianEigens(Configuration c, Potential * pPot, std::vector<double> & results, std::vector< std::vector<double> > * pEigenVectors = nullptr);

//special version for pair potential, and produce only lowest n eigen values for better efficiency
//not robust for very small systems (N<100), not sure why.
void LowestNHessianEigens(const Configuration & c, PairPotential & pot, size_t n, std::vector<double> & result2, std::vector< std::vector<double> > * pEigenVectors = nullptr);

#endif
