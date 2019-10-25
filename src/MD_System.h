#ifndef SYSTEM_INCLUDED
#define SYSTEM_INCLUDED

#include <functional>

#include "Potential.h"
#include "RandomGenerator.h"
#include "PeriodicCellList.h"




class ParticleMolecularDynamics
{
private:
	//before using, potential must be set to use the current configuration.
	//get Acceleration of all particles
	void _GetAcc(Potential & potential);
public:
	double TimeStep;
	double Time;//A Time Counter
	std::vector<double> Mass;
	Configuration Position;
	std::vector<GeometryVector> Velocity , Acceleration;
	size_t NumParticle;
	unsigned short Dimension;
	double MaxDeltaForce;//maximum delta F, user can query this variable to determine step size

	//variable for Nose-Hoover algorithm
	double xi;
	double MassQ;

	////////////////////////////////////////////////////////
	void InitNoseHoover(double Q);
	ParticleMolecularDynamics(const Configuration & Position, double TimeStep, double Mass);
	~ParticleMolecularDynamics();
	ParticleMolecularDynamics(std::istream & ifile);
	void WriteBinary(std::ostream & ofile);

	double GetPotentialEnergy(Potential & potential);
	void Evolve(size_t repeat, Potential & potential);
	void Evolve_AutoTimeStep(size_t repeat, Potential & potential, double DeltaLogEnergyPerStep);
	void Evolve_RecordDeltaForce(size_t repeat, Potential & potential);
	void NoseHooverEvolve(size_t repeat, Potential & potential, double kT);

	void AndersonEvolve(size_t repeat, Potential & potential, double kT, double ResetSpeedProbability, RandomGenerator & gen);
	double GetKineticEnergy(void);
	void SetRandomSpeed(double kT, RandomGenerator & gen);
	void RemoveOverallMomentum(void);

	//search for the optimal time step size
	//variable TimeStep as well as this->TimeStep was adjusted so that (maximum delta force)/(average force)
	//between two steps is close to TargetDeltaForce
	//
	//This function DESTROYS equilibrium. Re-equilibrate after calling this function
	void FindOptimalTimeStep(Potential & potential, double & TimeStep, double TargetDeltaForce, double kT);
};

// Basis vector moves are not implemented. i.e. simulation box is not changed
class ParticleMDAnnealing
{
public:
	ParticleMolecularDynamics * pSys;
	double alpha;
	double Tcut;
	double Tinit;
	double Pressure;
	RandomGenerator gen;
	Configuration * pList;
	std::ostream * Output;
	static bool StopSignal;
	bool SameLocalMinimaDetection;
	bool BoxMCMove;
	bool AutoTimeStep;
	bool Continue;//if true, read the last configuration from MDAnneal.ConfigPack and continue from there
					//also if true, will append to MDAnneal.ConfigPack rather than overwrite it.
	bool RelaxStructureBeforeAnnealing;
	std::string ConfigPackName;//default: "MDAnneal"

	std::function<void(Configuration & List, Potential & pot, double Pressure, double MinDistance)> RelaxStructureFunction;

	//ParticleMDAnnealing(DimensionType Dimension, size_t NumParticle, double InitialDensity, int RandomSeed, double Pressure, Potential & pot, std::ostream * output = nullptr, double CoolingScheduleRescale = 1, double MinDistance=0, bool UseSortedList=false);//AlgorithmLevel : 0:Structure Optimization, 1:Monte Carlo, >1:Parallel Tempering with n copies
	ParticleMDAnnealing(const Configuration & pInitConfig, int RandomSeed, double Pressure, Potential & pot, bool BoxMCMove = false, std::ostream * output = nullptr, double CoolingScheduleRescale = 1, bool UseSortedList=false);
	void Anneal(Potential & pot, std::function<void (void)> Callback=[](){});

	~ParticleMDAnnealing()
	{
		//if(this->pList!=nullptr)
		//	delete this->pList;
		if(this->pSys!=nullptr)
			delete this->pSys;
	}
};

#endif
