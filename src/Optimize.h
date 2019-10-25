#ifndef OPTIMIZE_INCLUDED
#define OPTIMIZE_INCLUDED

#include "LatticeSumSolver.h"
#include "ParameterSet.h"
#include <vector>
#include <stdexcept>
#include <string>

double GetPressure(LatticeSumSolver * l, double volume, IsotropicPairPotential & pot);//Get the equilibrium pressure of this structure at this specific volume
bool ReOptimize(ParameterSet * & pInputParam, const std::vector<LatticeSumSolver *> & Solvers, size_t ReOptimizeTime, double TargetStructureVolume, double & result, size_t & NumOutput, double & CenterPressure, double & DeltaPressure, RandomGenerator & gen, std::string OutputPrefix="");
double LowestEigenValueOfStiffness(Configuration list, Potential & pot, double Pressure);
double LowestEigenValueOfDeformations(Configuration list, Potential & pot, double Pressure);


class MinDistanceTooBig : public std::exception
{
public:
	MinDistanceTooBig()
	{
	}
	virtual const char * what() const throw()
	{
		return "MinDistance is too big, a Competitor have pair distance smaller than this!\n";
	}
};

namespace EnthalpyDifference
{
	//return the highest of (E_Target-E_Competitor)
	//if pLowestCompetitor is not nullptr, then the list it's pointing to will be changed to the lowest Competitor
	double EnthalpyDifference(ParameterSet & param, const std::vector<LatticeSumSolver *> & Solver, double TargetStructureVolume, double Pressure, Configuration * pLowestCompetitor = nullptr);

	//return true if the target is locally stable (relax to itself)
	bool LocallyStable(Configuration tar, Potential & Pot, double Pressure, double MinDistance);


	//SpecificStage==-1 : perform all stages of optimization
	//SpecificStage!=-1 : perform only one stage of optimization
	double Optimize(const std::vector<LatticeSumSolver *> & vpSolvers, ParameterSet & param, double & Pressure, double TargetVolume, int SpecificStage=-1);
};

void AllDerivatives(double *result, Potential * pPot, Configuration & tar, double Pressure);


#endif
