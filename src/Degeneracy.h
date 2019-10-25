#ifndef DEGENERACY_INCLUDED
#define DEGENERACY_INCLUDED

#include "LatticeSumSolver.h"
#include "RandomGenerator.h"
#include <vector>
#include "PeriodicCellList.h"
#include <cstring>

//find a structure with the same g2(r) as pReference->GetStructure()
void SearchDegenerate(LatticeSumSolver * pReference, RandomGenerator & gen, size_t SearchTime, std::vector<Configuration> & results, size_t NumParticle);

//find two structures with the same g2(r)
void SearchDegenerate2(RandomGenerator & gen, size_t SearchTime, std::vector<std::pair<Configuration, Configuration> > & results, size_t NumParticle1, size_t NumParticle2, DimensionType dim, double * pd2=nullptr, double * pd3=nullptr);

//find a structure with smallest number of particles, but has the same g2 and g3 than Reference
Configuration FindStructure(Configuration Reference, RandomGenerator & gen, size_t SearchTime=50, double Precision=1e-5);
Configuration FindStructure_SpecificNumberBasis(Configuration Reference, RandomGenerator & gen, size_t NumBasis, size_t SearchTime=50, double Precision=1e-5);

#endif