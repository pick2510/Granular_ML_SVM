#ifndef RandomSequentialAddition_Included
#define RandomSequentialAddition_Included 

#include <vector>
#include "PeriodicCellList.h"
#include <cstdint>
class VoxelList
{
private:
	typedef unsigned short voxeltype;
	int dimension, level;
	double halfsize, originhalfsize;
	std::vector<voxeltype> * originindex;
	std::vector<unsigned short> * subindex;
	unsigned short nbrsplitvoxel;
	signed char * splitvoxellist;

	//convert double * to GeometryVector
	GeometryVector Coord2Vect(const double * coord);
	void getsplitvoxellist(void);

public:
	VoxelList();
	~VoxelList();
	VoxelList(int Dimension, voxeltype VoxelRank, size_t ExpectedSize=10000);
	void GetOriginalList(const SpherePacking & Config, double Radius);
	unsigned long NumVoxel(void) const;
	void GetCenter(double * result, size_t nbr) const;//get the center coordinate of the nbr-th voxel
	void GetRandomCoordinate(double * result, RandomGenerator & gen) const;//generate a random coordinate which is inside a random voxel, write it into result
	bool CheckOverlapDeep(const SpherePacking & Config, const double & Radius, const double * Center, double QuaterSize, unsigned int checklevel);
	void SplitVoxelSerial(const SpherePacking & Config, double Radius, int checklevel=0);
	void DistillVoxelList(const SpherePacking & Config, double Radius, int checklevel=0);
	void SplitVoxel(const SpherePacking & Config, double Radius, int checklevel=0);
	double GetVoxelSize() const;
};

class RSA
{
public:
	SpherePacking * pSpheres;
	VoxelList * pVoxels;
	int dimension;
	RandomGenerator gen;
	//double diameter;
	double radius;
	double spherevolume;//the volume of a sphere

	//convert double * to GeometryVector
	GeometryVector Coord2Vect(const double * coord);

	//initialize square box with side length 1, sphere radii are set so that around NumSphere spheres exist in saturation limit
	RSA(int Dimension, unsigned long NumSphere);

	//initial configuration pa, upcoming spheres have one radii
	//pa should support spheres of different sizes, but this is not tested.
	RSA(const SpherePacking & pa, double radius);
	~RSA();

	void RSA_I(size_t NumInsertedLimit, size_t TrialTime);//try to insert sphere, if less than NumInsertedLimit spheres inserted in TrialTime trials, stop
	void RSA_I(size_t NumCut);//try to insert sphere until there are NumCut spheres
	void RSA_II_Serial( size_t NumInsertedLimit, size_t TrialTime);
	void RSA_II_Parallel( size_t NumInsertedLimit, size_t TrialTime);
	void RSA_II(size_t NumInsertedLimit, size_t TrialTime, bool parallel);
	double Density(void) const;
	void GetVoxelList(double VoxelDiameterRatio);//get original voxel list
	void SplitVoxel(bool parallel, int checklevel=0);
	void OutputSpheres(std::ostream & out);
};

//return the unsaturated packing with NumCut particles if it's achieved
SpherePacking GenerateRSAPacking(int dimension, unsigned long nbrsphere, short nbrinsertedlimit, unsigned long trial1, unsigned long trial2, double voxelratio, int seed, std::ostream & output, double * volumeratio=nullptr, bool Verbose=false, size_t NumCut=SIZE_MAX);


//return the unsaturated packing with NumCut particles if it's achieved
inline SpherePacking GenerateRSAPacking(int dimension, unsigned long nbrsphere, bool parallel=false, size_t NumCut=SIZE_MAX)
{
	unsigned long trial = 3*std::pow(10, dimension);
	return GenerateRSAPacking((-1)*std::abs(dimension), nbrsphere, parallel?(-3):3, trial, trial, 0.49, GetPreciseClock(), Verbosity>5?std::cout:logfile, nullptr, false, NumCut);
}
#endif