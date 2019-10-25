#ifndef MCSYSTEM_INCLUDED
#define MCSYSTEM_INCLUDED

#include <functional>
#include "MonteCarlo.h"
#include "PeriodicCellList.h"
#include "Move.h"


//construct this class, and set every data member if you want, then call Anneal().
//after that, get the configuration by COPYING (*pList)
class ParticleSimulatedAnnealing
{
private:
	MonteCarlo<Configuration, Potential *> * pSys;

public:
	double alpha;
	double Tcut;
	double Tinit;
	double CellMoveProbability;
	double Pressure;
	RandomGenerator gen;
	Configuration * pList;
	std::ostream * Output;
	double MinDistance;//abandon this move immediately if particles are within this distance
	static bool StopSignal;
	std::function<void(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance)> RelaxStructureFunction;

	ParticleSimulatedAnnealing(DimensionType Dimension, size_t NumParticle, double InitialDensity, int RandomSeed, double Pressure, Potential & pot, std::ostream * output = nullptr, double CoolingScheduleRescale = 1, double MinDistance=0, int AlgorithmLevel = 1, bool UseSortedList=false);//AlgorithmLevel : 0:Structure Optimization, 1:Monte Carlo, >1:Parallel Tempering with n copies
	ParticleSimulatedAnnealing(const Configuration & pInitConfig, int RandomSeed, double Pressure, Potential & pot, std::ostream * output = nullptr, double CoolingScheduleRescale = 1, double MinDistance=0, int AlgorithmLevel = 1, bool UseSortedList=false);
	void Anneal(Potential & pot);

	~ParticleSimulatedAnnealing()
	{
		if(this->pList!=nullptr)
			delete this->pList;
		if(this->pSys!=nullptr)
			delete this->pSys;
	}
};



#endif
