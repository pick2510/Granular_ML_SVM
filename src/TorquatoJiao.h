#ifndef TORQUATOJIAO_INCLUDED
#define TORQUATOJIAO_INCLUDED

//Anything related to linear programming application on jamming hard-sphere packings

#include "etc.h"
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include "PeriodicCellList.h"
#include <exception>


bool IsRattler_Count(DimensionType dim, const std::vector<GeometryVector> & contacts);
#ifdef USE_TJ
class LPClass;
class TJJammer
{
public:
	class error : public std::exception
	{
	public:
		std::string message;
		error(const char * a): message(a){}
		virtual const char* what() const throw()
		{
			return message.c_str();
		}
	};
	// input file variables - these are the defaults in case input isn't read!
	// See file "HopkinsTJ.standard.input" for details about all of these parameters
	std::string parentDirectory;
	std::string outputFilename;
	std::string tempFolder;
	bool   useNNL ;
	bool autoSave  ;//Whether to save progress in a temporary folder periodically
	int    N  ;                  //The number of spheres
	int    dim  ;                //The dimension of the problem
	bool   maxInitRad  ;     //If true, the radii of the spheres are maximized in the IC (so that the first contact is made)
	double transMax  ;     //The maximum translation of a sphere in a single LP step (lattice vectors?)
	double compMax  ;      //The maximum bulk (compressive) strain in a single step
	double shearMax  ;    //The maximum shear              strain in a single step
	int    strictJam  ;       //1 allows for full box deformations (as allowed by compMax and shearMax), but 0 resticts the compression to be shape-preserving (which will cause collective jamming instead of strict jamming)
	double delta  ;         //The influence sphere radius (?)
	int    maxIters  ;     //The maximum number of times the LP solver may be invoked
	std::string termCriterion;        //Either "maxDisp" or "latticeVol"
	double phiTerm  ;       //The maximum density allowed by the simulation (nice when we want to stop early)
	double termTol  ;    //The termination criterion: if the density changes by less than this, then we are done (?)
	//Alternatively, if the biggest sphere movement (relative to the center of gravity) was smaller than this, then we're done
	double MCStartSweeps  ; //How many MC sweeps to do before the SLP algorithm
	double MCMidSweeps  ;   //How many MC sweeps to do after each LP iteration          
	double mc_dispLimit;         //The displacement limit for a trial move in the MC routine
	double maxDisp;              //The farthest (relative to the total packing displacement) that a sphere has moved in an iteration
	bool   randLatVecs  ;   //When generating a packing, whether or not the lattice vectors should be random (otherwise, cube)
	double maxLatLength  ;  //The maximum length of a lattice vector (at initialization) (?)
	double minLatLength  ; //The minimum length of a lattice vector (at initialization)
	double minAngDeg  ;    //The minimum angle between vectors when starting the lattice.
	double manualMinVol  ; //
	double radLatTol ;
	int    topFail ;
	int    maxVolNum ;
	double inputPhi  ;        //The intial (starting) density
	int    overBoxes  ;          //How many rings (of periodic image cells) to check beyond the actual cell (0 = use L/2 method)
	int    maxBoxNum ;
	bool   randOverlapMove ;
	double resizeSpace ;
	double resizeTol ;
	double NNLextraDist ;
	int    printEvery  ;        //if >0, Write out the configuration to file after every [this many] iterations
	bool   printThisIter_pkg;
	bool   printThisIter_LP;
	int    printPrecision  ;   //Number of decimal places - outFile.precision(printPrecision);
	bool   useTimer  ;          //Whether to time the run (1) or not (0).  Units are microseconds
	bool   scaleUnitRadius  ;   //If yes, scale the packing so that the smalest sphere has radius = 1
	double rescaleIC  ;       //After you make/load an IC, multiply the size of the packing by this (to "rewind" packings)
	int    printLPEvery  ;      //if >0, Write out the LP's solution to file after every [this many] iterations


	//Solver variable:
	std::string solver  ;                //"GLPK" or "Gurobi"

	//GLPK solver params
	double feasibleTol ;
	double pivotTol  ;//Default value used by GLPK, 0<pt<1...increase to improve numerical accuracy...
	int GLPK_PresolveVar ;
	int GLPK_BasisFactType ;

	//Gurobi solver params
	std::string GRB_Method  ;//Which solver to use
	int    GRB_Presolve  ;   //How much presolving to do (-1,0,1,2) = (auto,none,standard,aggresive)
	int    grbThreads ;      //Limit how many threads Gurobi can use at once

	TJJammer();

protected:
	// Timer variables:
	class TimerClass 
	{
	public:
		double timePassed;

		void Start();
		void End();
		void Display();
	private:
		time_t startTime;
		time_t endTime;
	};
	class NeighborNode 
	{
	public:
		NeighborNode *next;
		int index;

		~NeighborNode() 
		{
			if (next != 0) 
			{
				delete next;
			}
		}
	};
	// input file variables - these are the defaults in case input isn't read!
	const gsl_rng_type * RNGtype;
	gsl_rng *RNG;
	TimerClass timer; //See timer.h for structure
	TimerClass saveTimer; //To save the data periodically

	// running variables
	gsl_matrix *lambdas;
	gsl_matrix *inverseLambdas;
	double *localCoords ;
	double *globalCoords ;
	double *radii ;
	double *distTempL ;
	double *distTempG ;
	NeighborNode **neighborListDelta ;
	NeighborNode **neighborListOverlap ;
	int lpIters  ;//The number of LP iterations that have gone by
	double maxSingleMove ;
	double biggestRad ;
	double getGlobalLength(double *inputVec) ;
	double getVol();
	int resizeIfOverlap(bool randOver) ;
	double getSphereVol();
	void RescalePacking();
	void DeleteOldNNL();
	int calcNNLs() ;
	double InternalTransMax  ;     //The maximum translation of a sphere in a single LP step (lattice vectors?)
	double InternalCompMax  ;      //The maximum bulk (compressive) strain in a single step
	double InternalShearMax  ;    //The maximum shear              strain in a single step
	double InternalDelta;
	double InternalNNLextraDist;
	void MaximizeRadii();
	void MC_AdjustDispLimit(double accRate);
	double MC_Move(int i);
	double MC_Sweep(int numTries);
	void MC_Init();
	int printPacking(std::string thisOutputFilename, bool addIterNum);
	bool CheckTermination(double lastLastVol, double latVol, double phi);
	void PrintTime(TimerClass timer);
	void ClearSaveData();
	void SaveProgress();
	SpherePacking result;
public:
	void TJIterations(void);//run iterations until stop criteria is met
	void SetPacking(const SpherePacking & p);
	SpherePacking GetPacking(void);
	friend class LPClass;
};

//find if a sphere is a ratteler or not given it has these contacts to non-ratteler spheres
//need to loop through all spheres until no extra rattler are found
//This function used GLPK, which is NOT thread-safe
bool IsRattler_LP(DimensionType dim, const std::vector<GeometryVector> & contacts);

#else
class TJJammer
{
public:
	// input file variables - these are the defaults in case input isn't read!
	// See file "HopkinsTJ.standard.input" for details about all of these parameters
	std::string parentDirectory;
	std::string outputFilename;
	std::string tempFolder;
	bool   useNNL ;
	bool autoSave  ;//Whether to save progress in a temporary folder periodically
	int    N  ;                  //The number of spheres
	int    dim  ;                //The dimension of the problem
	bool   maxInitRad  ;     //If true, the radii of the spheres are maximized in the IC (so that the first contact is made)
	double transMax  ;     //The maximum translation of a sphere in a single LP step (lattice vectors?)
	double compMax  ;      //The maximum bulk (compressive) strain in a single step
	double shearMax  ;    //The maximum shear              strain in a single step
	int    strictJam  ;       //1 allows for full box deformations (as allowed by compMax and shearMax), but 0 resticts the compression to be shape-preserving (which will cause collective jamming instead of strict jamming)
	double delta  ;         //The influence sphere radius (?)
	int    maxIters  ;     //The maximum number of times the LP solver may be invoked
	std::string termCriterion;        //Either "maxDisp" or "latticeVol"
	double phiTerm  ;       //The maximum density allowed by the simulation (nice when we want to stop early)
	double termTol  ;    //The termination criterion: if the density changes by less than this, then we are done (?)
	//Alternatively, if the biggest sphere movement (relative to the center of gravity) was smaller than this, then we're done
	double MCStartSweeps  ; //How many MC sweeps to do before the SLP algorithm
	double MCMidSweeps  ;   //How many MC sweeps to do after each LP iteration          
	double mc_dispLimit;         //The displacement limit for a trial move in the MC routine
	double maxDisp;              //The farthest (relative to the total packing displacement) that a sphere has moved in an iteration
	bool   randLatVecs  ;   //When generating a packing, whether or not the lattice vectors should be random (otherwise, cube)
	double maxLatLength  ;  //The maximum length of a lattice vector (at initialization) (?)
	double minLatLength  ; //The minimum length of a lattice vector (at initialization)
	double minAngDeg  ;    //The minimum angle between vectors when starting the lattice.
	double manualMinVol  ; //
	double radLatTol ;
	int    topFail ;
	int    maxVolNum ;
	double inputPhi  ;        //The intial (starting) density
	int    overBoxes  ;          //How many rings (of periodic image cells) to check beyond the actual cell (0 = use L/2 method)
	int    maxBoxNum ;
	bool   randOverlapMove ;
	double resizeSpace ;
	double resizeTol ;
	double NNLextraDist ;
	int    printEvery  ;        //if >0, Write out the configuration to file after every [this many] iterations
	bool   printThisIter_pkg;
	bool   printThisIter_LP;
	int    printPrecision  ;   //Number of decimal places - outFile.precision(printPrecision);
	bool   useTimer  ;          //Whether to time the run (1) or not (0).  Units are microseconds
	bool   scaleUnitRadius  ;   //If yes, scale the packing so that the smalest sphere has radius = 1
	double rescaleIC  ;       //After you make/load an IC, multiply the size of the packing by this (to "rewind" packings)
	int    printLPEvery  ;      //if >0, Write out the LP's solution to file after every [this many] iterations


	//Solver variable:
	std::string solver  ;                //"GLPK" or "Gurobi"

	//GLPK solver params
	double feasibleTol ;
	double pivotTol  ;//Default value used by GLPK, 0<pt<1...increase to improve numerical accuracy...
	int GLPK_PresolveVar ;
	int GLPK_BasisFactType ;

	//Gurobi solver params
	std::string GRB_Method  ;//Which solver to use
	int    GRB_Presolve  ;   //How much presolving to do (-1,0,1,2) = (auto,none,standard,aggresive)
	int    grbThreads ;      //Limit how many threads Gurobi can use at once

	TJJammer()
	{
		std::cerr<<"Error in TJJammer : TJ is not enabled!\n";
	}

protected:
public:
	void TJIterations(void)//run iterations until stop criteria is met
	{
		std::cerr<<"Error in TJJammer : TJ is not enabled!\n";
	}
	void SetPacking(const SpherePacking & p)
	{
		std::cerr<<"Error in TJJammer : TJ is not enabled!\n";
	}
	SpherePacking GetPacking(void)
	{
		std::cerr<<"Error in TJJammer : TJ is not enabled!\n";
		return SpherePacking();
	}
};
inline bool IsRattler_LP(DimensionType dim, const std::vector<GeometryVector> & contacts)
{
	std::cerr << "Error in IsRattler_LP : TJ is not enabled!\n"; 
		return false;
}
#endif

#endif