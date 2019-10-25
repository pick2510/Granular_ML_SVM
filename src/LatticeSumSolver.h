#ifndef LATTICESUMSOLVER_INCLUDED
#define LATTICESUMSOLVER_INCLUDED

#include <vector>
#include "Potential.h"
#include "PeriodicCellList.h"

struct LatticeTerm
{
	unsigned long number;
	double distance;
};



//partially deprecated class

//Initially, this class was used to perform lattice sum. the class can generate the theta series to a given distance
//and calculate the lattice sum
//However, when we start supporting multicomponent system, this class has not been rewritten to include multiple theta series.

//Therefore, the member function LatticeSum() is deprecate. To revive the class, rewrite it to include multiple theta series for each pair of component.
//This class is only used for SolverDistance().

class LatticeSumSolver
{
public:
	std::vector<LatticeTerm> Terms;//Terms per unit cell
	double CurrentRc;
	double OriginalDensity;//the density of the crystal that generates Terms
	size_t NumParticle;//Number of atoms per unit cell. 
	DimensionType Dimension;

	LatticeSumSolver(void);
	virtual const char * Tag(void) =0;
	virtual ~LatticeSumSolver();
	virtual Configuration GetStructure(void) =0;
	virtual void UpdateTerms(double NewRc);
	double LatticeSum(double Volume, IsotropicPairPotential & potential);
};

bool SameSolver(LatticeSumSolver * a, LatticeSumSolver * b);
double SolverDistance(LatticeSumSolver * a, LatticeSumSolver * b);

#endif