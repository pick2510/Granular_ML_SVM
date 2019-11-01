#include <fstream>
#include <iostream>
#include <complex>
#include <list>
#include <set>
#include <numeric>
#include <memory>
#include <future>
#include <boost/math/special_functions.hpp>
#ifdef USE_EIGEN
#include <Eigen/Eigen>
#endif
#include "PeriodicCellList.h"
#include "Plots.h"
#include "MC_System.h"
#include "MD_System.h"
#include "StructureOptimization.h"
#include "PhononFrequency.h"
#include "Voronoi.h"
#include "Optimize.h"
#include "Solvers.h"
#include "ThreeBody.h"
#include "Degeneracy.h"
#include "TorquatoJiao.h"

#ifdef USE_EIGEN
//calculate d2min between configurations now and next for particle i
template <int Dimension>
double D2Min_Inner(const Configuration &now, const Configuration &next, size_t i, double Cutoff, double *pAffinePart = nullptr, double *pShearStrain = nullptr, double *pIsotropicStrain = nullptr)
{
	//Configuration is a box of arbitrary shape and particles with a cell list implemented
	double DispSqSum;

	Eigen::Matrix<double, Dimension, Dimension> X;
	Eigen::Matrix<double, Dimension, Dimension> Y;
	Eigen::Matrix<double, Dimension, Dimension> Transform;

	double D2Min = 0.0;
	double Neighbors = 0.0;
	X = Eigen::Array<double, Dimension, Dimension>::Zero();
	Y = Eigen::Array<double, Dimension, Dimension>::Zero();
	DispSqSum = 0.0;
	Transform = Eigen::Array<double, Dimension, Dimension>::Zero();

	//GeometryVector is a d-dimensional vector
	auto IterateFunc = [&](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t Sourceparticle) {
		if (Sourceparticle != i)
		{
			//Relative coordinates are coordinates before applying CBox::Transform()
			GeometryVector dj = next.GetRelativeCoordinates(Sourceparticle) - next.GetRelativeCoordinates(i);
			next.RelativeCoordToMinimumImage(dj);
			//convert to Cartesian coordinates, which are coordinates after Transform
			dj = next.RelativeCoord2CartesianCoord(dj);

			Eigen::Matrix<double, Dimension, 1> displacement_j(dj.x), displacement(shift.x);

			X += displacement_j * displacement.transpose();
			Y += displacement * displacement.transpose();
			DispSqSum += dj.Modulus2();

			Neighbors++;

			//debug temp
			//std::cout <<Sourceparticle << "\t" << shift;
		}
	};
	//Iterate through neighbors of particle i, for each neighbor, call IterateFunc()
	now.IterateThroughNeighbors(i, Cutoff, IterateFunc);
	//debug temp
	//std::cout << "Neighbors=" << Neighbors << std::endl;

	Transform = X * Y.inverse();
	double affinePart = (-2.0 * (Transform * X.transpose()).trace() + (Transform * Y * Transform.transpose()).trace()) * (-1.0);
	if (pAffinePart != nullptr)
		(*pAffinePart) = affinePart;
	if (pIsotropicStrain != nullptr)
	{
		double ave = Transform.trace() / Dimension;
		(*pIsotropicStrain) = ave - 1.0;
	}
	if (pShearStrain != nullptr)
	{
		for (int i = 0; i + 1 < Dimension; i++)
			for (int j = i + 1; j < Dimension; j++)
			{
				double ave = 0.5 * (Transform(i, j) + Transform(j, i));
				Transform(i, j) = ave;
				Transform(j, i) = ave;
			}
		Eigen::SelfAdjointEigenSolver<decltype(Transform)> solver(Transform);
		if (solver.info() != Eigen::Success)
		{
			std::cerr << "Unsuccessful eigenvalue calculation\n";
			(*pShearStrain) = 0.0;
		}
		else
		{
			auto evalues = solver.eigenvalues();
			double sum = 0.0;
			for (int i = 0; i < Dimension; i++)
				sum += evalues(i);
			double ave = sum / Dimension;
			if (pShearStrain != nullptr)
			{
				(*pShearStrain) = 0.0;
				for (int i = 0; i < Dimension; i++)
					(*pShearStrain) += std::abs(evalues(i) - ave);
			}
		}
	}
	return (DispSqSum - affinePart) / Neighbors;
}
#else
template <int Dimension>
double D2Min_Inner(const Configuration &now, const Configuration &next, size_t i, double Cutoff, double *pAffinePart = nullptr, double *pShearStrain = nullptr, double *pIsotropicStrain = nullptr)
{
	std::cerr << "Error in D2Min_Inner : Eigen is not enabled\n";
	return 0.0;
}
#endif

double D2Min(const Configuration &now, const Configuration &next, size_t i, double Cutoff, double *pAffinePart = nullptr, double *pShearStrain = nullptr, double *pIsotropicStrain = nullptr)
{
	switch (now.GetDimension())
	{
	case 1:
		return D2Min_Inner<1>(now, next, i, Cutoff, pAffinePart, pShearStrain, pIsotropicStrain);
	case 2:
		return D2Min_Inner<2>(now, next, i, Cutoff, pAffinePart, pShearStrain, pIsotropicStrain);
	case 3:
		return D2Min_Inner<3>(now, next, i, Cutoff, pAffinePart, pShearStrain, pIsotropicStrain);
	case 4:
		return D2Min_Inner<4>(now, next, i, Cutoff, pAffinePart, pShearStrain, pIsotropicStrain);
	default:
		std::cerr << "Error in D2Min : unsupported dimension!\n";
		return 0.0;
	}
}

#include "DisplaySpheres.h"
//several color bars
//gk: green-black
//rycbk: red-yellow-cyan-blue-black
//Rrycb: dark red-red-yellow-cyan-blue
void AssignColor_gk(double ColorIndex, sphere &temp)
{
	if (ColorIndex > 1.0)
	{
		temp.red = 0.0;
		temp.green = 1.0;
		temp.blue = 0.0;
	}
	else if (ColorIndex > 0.0)
	{
		double d = (ColorIndex - 0.0) / 1.0;
		temp.red = 0.0;
		temp.green = d;
		temp.blue = 0.0;
	}
	else
	{
		temp.red = 0.0;
		temp.green = 0.0;
		temp.blue = 0.0;
	}
	temp.transparency = 1.00;
}
void AssignColor_rycbk(double ColorIndex, sphere &temp)
{
	if (ColorIndex > 1.0)
	{
		temp.red = 1.0;
		temp.green = 0.0;
		temp.blue = 0.0;
	}
	else if (ColorIndex > 0.75)
	{
		double d = (ColorIndex - 0.75) / 0.25;
		temp.red = 1.0;
		temp.green = 1.0 - d;
		temp.blue = 0.0;
	}
	else if (ColorIndex > 0.5)
	{
		double d = (ColorIndex - 0.5) / 0.25;
		temp.red = d;
		temp.green = 1.0;
		temp.blue = 1.0 - d;
	}
	else if (ColorIndex > 0.25)
	{
		double d = (ColorIndex - 0.25) / 0.25;
		temp.red = 0.0;
		temp.green = d;
		temp.blue = 1.0;
	}
	else if (ColorIndex > 0.0)
	{
		double d = (ColorIndex - 0.0) / 0.25;
		temp.red = 0.0;
		temp.green = 0.0;
		temp.blue = d;
	}
	else
	{
		temp.red = 0.0;
		temp.green = 0.0;
		temp.blue = 0.0;
	}
	temp.transparency = 1.00;
}
void AssignColor_Rrycb(double ColorIndex, sphere &temp)
{
	if (ColorIndex > 1.0)
	{
		temp.red = 0.5;
		temp.green = 0.0;
		temp.blue = 0.0;
	}
	else if (ColorIndex > 0.75)
	{
		double d = (ColorIndex - 0.75) / 0.25;
		temp.red = 1.0 - 0.5 * d;
		temp.green = 0.0;
		temp.blue = 0.0;
	}
	else if (ColorIndex > 0.5)
	{
		double d = (ColorIndex - 0.5) / 0.25;
		temp.red = 1.0;
		temp.green = 1.0 - d;
		temp.blue = 0.0;
	}
	else if (ColorIndex > 0.25)
	{
		double d = (ColorIndex - 0.25) / 0.25;
		temp.red = d;
		temp.green = 1.0;
		temp.blue = 1.0 - d;
	}
	else if (ColorIndex > 0.0)
	{
		double d = (ColorIndex - 0.0) / 0.25;
		temp.red = 0.0;
		temp.green = d;
		temp.blue = 1.0;
	}
	else
	{
		temp.red = 0.0;
		temp.green = 0.0;
		temp.blue = 1.0;
	}
	temp.transparency = 1.00;
}

void AssignColor_kw(double ColorIndex, sphere &temp)
{
	if (ColorIndex > 1.0)
	{
		temp.red = 0.0;
	}
	else if (ColorIndex > 0)
	{
		temp.red = 1.0 - ColorIndex;
	}
	else
	{
		temp.red = 1.0;
	}
	temp.transparency = 1.0;
	temp.transparency = 1.0 - temp.red;
	temp.red = 0.0;
	temp.green = temp.red;
	temp.blue = temp.red;
}

void DisplaySSPackingMovie_SingleColorIndex(ConfigurationPack &pk, double RadiusA, double RadiusB, std::function<double(long, long)> GetColorIndexFunc, std::function<void(double, sphere &)> AssignColorFunc, long numFrames)
{
	if (pk.NumConfig() == 0)
	{
		std::cout << "warning in Display3DConfigurationMovie : ConfigurationPack is empty. Returning!\n";
		return;
	}
	Configuration a = pk.GetConfig(0);

	auto GetFrameFunc = [&](std::vector<sphere> &spheres, std::vector<line> &lines, long FrameNumber) -> void {
		spheres.clear();
		lines.clear();
		Configuration c1 = pk.GetConfig(FrameNumber);

		ConfigurationToSpheres(c1, spheres, lines, RadiusA);

		for (size_t i = 0; i < c1.NumParticle(); i++)
		{
			sphere &temp = spheres[i];

			char chr1 = c1.GetCharacteristics(i).name[0];

			if (chr1 != 'A')
				temp.radius *= RadiusB / RadiusA;

			double ColorIndex = GetColorIndexFunc(FrameNumber, i);
			AssignColorFunc(ColorIndex, temp);
		}
	};

	DisplaySpheres_Movie(GetFrameFunc, numFrames);
}

//display a 2D configuration by a green ring and a magenta core, controlled by GetColor1Func and GetColor2Func
void DisplaySSPackingMovie_DoubleColorIndex_green_magenta(ConfigurationPack &pk, double RadiusA, double RadiusB, std::function<double(long, long)> GetColor1Func, std::function<double(long, long)> GetColor2Func, long numFrames)
{
	::OrthoProjection = true;
	if (pk.NumConfig() == 0)
	{
		std::cout << "warning in Display3DConfigurationMovie : ConfigurationPack is empty. Returning!\n";
		return;
	}
	Configuration a = pk.GetConfig(0);

	auto GetFrameFunc = [&](std::vector<sphere> &spheres, std::vector<line> &lines, long FrameNumber) -> void {
		spheres.clear();
		lines.clear();
		Configuration c1 = pk.GetConfig(FrameNumber);

		ConfigurationToSpheres(c1, spheres, lines, RadiusA);

		for (size_t i = 0; i < c1.NumParticle(); i++)
		{
			sphere &temp = spheres[i];

			char chr1 = c1.GetCharacteristics(i).name[0];

			if (chr1 != 'A')
				temp.radius *= RadiusB / RadiusA;

			double Color1Index = GetColor1Func(FrameNumber, i);
			double Color2Index = GetColor2Func(FrameNumber, i);
			AssignColor_gk(Color1Index, spheres[i]);

			spheres.push_back(spheres[i]);
			AssignColor_gk(Color2Index, spheres.back());
			spheres.back().radius *= 0.5;
			spheres.back().z -= 0.1;
			std::swap(spheres.back().green, spheres.back().red);
			spheres.back().blue = spheres.back().red;
		}
	};
	DisplaySpheres_Movie(GetFrameFunc, numFrames);
}

void DisplaySSPackingMovie_DoubleColorIndex_blackCircleThickness_rainbow(ConfigurationPack &pk, double RadiusA, double RadiusB, std::function<double(long, long)> GetColor1Func, std::function<double(long, long)> GetColor2Func, long numFrames)
{
	::OrthoProjection = true;
	if (pk.NumConfig() == 0)
	{
		std::cout << "warning in Display3DConfigurationMovie : ConfigurationPack is empty. Returning!\n";
		return;
	}
	Configuration a = pk.GetConfig(0);

	auto GetFrameFunc = [&](std::vector<sphere> &spheres, std::vector<line> &lines, long FrameNumber) -> void {
		spheres.clear();
		lines.clear();
		Configuration c1 = pk.GetConfig(FrameNumber);

		ConfigurationToSpheres(c1, spheres, lines, RadiusA);

		for (size_t i = 0; i < c1.NumParticle(); i++)
		{
			sphere &temp = spheres[i];

			char chr1 = c1.GetCharacteristics(i).name[0];

			if (chr1 != 'A')
				temp.radius *= RadiusB / RadiusA;

			double Color1Index = GetColor1Func(FrameNumber, i);
			double Color2Index = GetColor2Func(FrameNumber, i);
			spheres[i].red = 0.0;
			spheres[i].green = 0.0;
			spheres[i].blue = 0.0;
			spheres[i].transparency = 1.0;

			spheres.push_back(spheres[i]);
			AssignColor_Rrycb(Color2Index, spheres.back());
			if (Color1Index > 1.0)
			{
				spheres[i].radius *= 2.0;
				spheres.back().radius *= 0.5;
			}
			else if (Color1Index > 0.0)
			{
				spheres[i].radius *= 1.0 + Color1Index;
				spheres.back().radius *= 1.0 - 0.5 * Color1Index;
			}
			spheres.back().z -= 0.1;
		}
	};
	DisplaySpheres_Movie(GetFrameFunc, numFrames);
}

void DisplaySSPackingMovie_DoubleColorIndex_transparency_rainbow(ConfigurationPack &pk, double RadiusA, double RadiusB, std::function<double(long, long)> GetColor1Func, std::function<double(long, long)> GetColor2Func, long numFrames)
{
	::OrthoProjection = true;
	if (pk.NumConfig() == 0)
	{
		std::cout << "warning in Display3DConfigurationMovie : ConfigurationPack is empty. Returning!\n";
		return;
	}
	Configuration a = pk.GetConfig(0);

	auto GetFrameFunc = [&](std::vector<sphere> &spheres, std::vector<line> &lines, long FrameNumber) -> void {
		spheres.clear();
		lines.clear();
		Configuration c1 = pk.GetConfig(FrameNumber);

		ConfigurationToSpheres(c1, spheres, lines, RadiusA);

		for (size_t i = 0; i < c1.NumParticle(); i++)
		{
			sphere &temp = spheres[i];

			char chr1 = c1.GetCharacteristics(i).name[0];

			if (chr1 != 'A')
				temp.radius *= RadiusB / RadiusA;

			double Color1Index = GetColor1Func(FrameNumber, i);
			double Color2Index = GetColor2Func(FrameNumber, i);
			AssignColor_Rrycb(Color2Index, spheres[i]);
			spheres[i].transparency = Color1Index;
		}
	};
	DisplaySpheres_Movie(GetFrameFunc, numFrames);
}

void DisplaySSPackingMovie_D2Min(ConfigurationPack &pk, double RadiusA, double RadiusB, double D2MinCutoff, double D2MinScale = 0.2, size_t D2min_Interval = 1)
{
	Configuration c0 = pk.GetConfig(0);
	//pre-compute d2min
	size_t nParticle = c0.NumParticle();
	std::vector<double> d2min(nParticle * (pk.NumConfig() - D2min_Interval), 0.0);
#pragma omp parallel for schedule(static)
	for (int f = 0; f < pk.NumConfig() - D2min_Interval; f++)
	{
		Configuration temp = pk.GetConfig(f);
		Configuration temp1 = pk.GetConfig(f + D2min_Interval);
		for (int n = 0; n < temp.NumParticle(); n++)
			d2min[f * nParticle + n] = D2Min(temp, temp1, n, D2MinCutoff) / D2MinScale;
	}
	DisplaySSPackingMovie_SingleColorIndex(pk, RadiusA, RadiusB, [&](size_t f, size_t n) -> double {
		return d2min[f * nParticle + n];
	},
										   AssignColor_rycbk, pk.NumConfig() - D2min_Interval);
}

void DisplaySSPackingMovie_D2Min(ConfigurationPack &pk, double RadiusA, double RadiusB, const std::vector<double> &d2min, double D2MinScale)
{
	size_t np = pk.GetConfig(0).NumParticle();
	DisplaySSPackingMovie_SingleColorIndex(pk, RadiusA, RadiusB, [&](size_t f, size_t n) -> double {
		return d2min[f * np + n] / D2MinScale;
	},
										   AssignColor_rycbk, d2min.size() / np);
}

template <typename T>
void DumpStdVector(const std::vector<T> &data, const std::string &filename)
{
	std::fstream ofile(filename, std::fstream::out | std::fstream::binary);
	long long len = data.size();
	ofile.write((char *)(&len), sizeof(long long));
	if (len > 0)
		ofile.write((char *)(data.data()), len * sizeof(T));
}
template <typename T>
std::vector<T> ReadStdVectorDump(const std::string &filename)
{
	std::fstream ifile(filename, std::fstream::in | std::fstream::binary);
	if (!ifile.good())
		return std::vector<T>();
	long long len;
	ifile.read((char *)(&len), sizeof(long long));
	if (ifile.eof())
		return std::vector<T>();
	std::vector<T> result;
	result.resize(len);
	ifile.read((char *)(result.data()), len * sizeof(T));
	return result;
}

void RemoveRattler(Configuration &c, PairPotential &Pot)
{
	std::vector<std::vector<GeometryVector>> allContacts;
	std::vector<std::vector<size_t>> allIndex;
	std::set<size_t> toConsider;
	std::set<size_t> toRemove;
	for (size_t i = 0; i < c.NumParticle(); i++)
	{
		std::vector<GeometryVector> contacts;
		std::vector<size_t> index;
		AtomInfo a1 = c.GetCharacteristics(i);
		c.IterateThroughNeighbors(i, Pot.Rcut, [&](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t Sourceparticle) -> void {
			if (Sourceparticle != i || LatticeShift.Modulus2() != 0.0)
			{
				AtomInfo a2 = c.GetCharacteristics(Sourceparticle);
				GeometryVector f;
				Pot.ForceFormula(shift, f, a1, a2);
				if (f.Modulus2() != 0.0)
				{
					contacts.push_back(shift);
					index.push_back(Sourceparticle);
				}
			}
		});
		allContacts.push_back(contacts);
		allIndex.push_back(index);
	}
	for (size_t i = 0; i < c.NumParticle(); i++)
		toConsider.insert(i);

	while (!toConsider.empty())
	{
		size_t i = *toConsider.begin();
		toConsider.erase(i);
		if (toRemove.find(i) == toRemove.end())
		{
			//not already removed, lets consider it
			//first, find out contacts that has not been removed
			for (size_t j = 0; j < allIndex[i].size(); j++)
			{
				if (toRemove.find(allIndex[i][j]) != toRemove.end())
				{
					allIndex[i][j] = allIndex[i].back();
					allIndex[i].pop_back();
					allContacts[i][j] = allContacts[i].back();
					allContacts[i].pop_back();
					j--;
				}
			}
			//second, LP test
			if (IsRattler_LP(c.GetDimension(), allContacts[i]))
			{
				toRemove.insert(i);
				for (size_t j = 0; j < allIndex[i].size(); j++)
					toConsider.insert(allIndex[i][j]);
			}
		}
	}
	for (auto iter = toRemove.rbegin(); iter != toRemove.rend(); ++iter)
	{
		c.DeleteParticle(*iter);
	}
}

//same as RemoveRattler, but whenever removing a particle in c, also remove the corresponding one in c1
void RemoveRattler2(Configuration &c, Configuration &c1, PairPotential &Pot)
{
	std::vector<std::vector<GeometryVector>> allContacts;
	std::vector<std::vector<size_t>> allIndex;
	std::set<size_t> toConsider;
	std::set<size_t> toRemove;
	for (size_t i = 0; i < c.NumParticle(); i++)
	{
		std::vector<GeometryVector> contacts;
		std::vector<size_t> index;
		AtomInfo a1 = c.GetCharacteristics(i);
		c.IterateThroughNeighbors(i, Pot.Rcut, [&](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t Sourceparticle) -> void {
			if (Sourceparticle != i || LatticeShift.Modulus2() != 0.0)
			{
				AtomInfo a2 = c.GetCharacteristics(Sourceparticle);
				GeometryVector f;
				Pot.ForceFormula(shift, f, a1, a2);
				if (f.Modulus2() != 0.0)
				{
					contacts.push_back(shift);
					index.push_back(Sourceparticle);
				}
			}
		});
		allContacts.push_back(contacts);
		allIndex.push_back(index);
	}
	for (size_t i = 0; i < c.NumParticle(); i++)
		toConsider.insert(i);

	while (!toConsider.empty())
	{
		size_t i = *toConsider.begin();
		toConsider.erase(i);
		if (toRemove.find(i) == toRemove.end())
		{
			//not already removed, lets consider it
			//first, find out contacts that has not been removed
			for (size_t j = 0; j < allIndex[i].size(); j++)
			{
				if (toRemove.find(allIndex[i][j]) != toRemove.end())
				{
					allIndex[i][j] = allIndex[i].back();
					allIndex[i].pop_back();
					allContacts[i][j] = allContacts[i].back();
					allContacts[i].pop_back();
					j--;
				}
			}
			//second, LP test
			if (IsRattler_LP(c.GetDimension(), allContacts[i]))
			{
				toRemove.insert(i);
				for (size_t j = 0; j < allIndex[i].size(); j++)
					toConsider.insert(allIndex[i][j]);
			}
		}
	}
	for (auto iter = toRemove.rbegin(); iter != toRemove.rend(); ++iter)
	{
		c.DeleteParticle(*iter);
		c1.DeleteParticle(*iter);
	}

	//std::cout << c.NumParticle()<< std::endl;
}

void DisplaySSPackingMovie_D2Min_softness(ConfigurationPack &pk, const std::vector<float> &softness, double RadiusA, double RadiusB, double D2MinScale, double D2MinCutoff, double SoftnessScale, size_t D2min_Interval)
{
	if (pk.NumConfig() <= D2min_Interval)
	{
		std::cout << "warning in Display3DConfigurationMovie : ConfigurationPack is too small. Returning!\n";
		return;
	}

	double meanS = 0.0;
	for (auto &s : softness)
		meanS += s;
	meanS /= softness.size();

	Configuration temp = pk.GetConfig(0);
	Configuration temp1 = pk.GetConfig(D2min_Interval);
	size_t currentf = 0;
	size_t softnessReadyF = pk.NumConfig();
	auto getD2minFunc = [&](size_t f, size_t n) -> double {
		if (f != currentf)
		{
			temp = pk.GetConfig(f);
			temp1 = pk.GetConfig(f + D2min_Interval);
			currentf = f;
		}
		return D2Min(temp, temp1, n, D2MinCutoff) / D2MinScale;
	};
	auto getSoftnessFunc = [&](size_t f, size_t n) -> double {
		return (softness[n + f * temp.NumParticle()] - meanS + SoftnessScale) / (2.0 * SoftnessScale);
	};
	DisplaySSPackingMovie_DoubleColorIndex_transparency_rainbow(pk, RadiusA, RadiusB, getD2minFunc, getSoftnessFunc, pk.NumConfig() - D2min_Interval);
}

Configuration readOmidPacking(std::fstream &ifile)
{
	if (!ifile.good())
		return Configuration();
	GeometryVector bas[3] = {GeometryVector(0.65, 0.0, 0.0), GeometryVector(0.0, 0.4, 0.0), GeometryVector(0.0, 0.0, 0.4)};
	Configuration pak(3, bas, 0.008, false);
	struct particle
	{
		int id;
		GeometryVector car;
		double radius;
	};
	std::vector<particle> ps;
	for (;;)
	{
		particle p;
		ifile >> p.id;

		double temp;
		ifile >> temp;

		p.car.SetDimension(3);
		ifile >> p.car.x[0] >> p.car.x[1] >> p.car.x[2];
		ifile >> p.radius;

		ps.push_back(p);

		ifile.ignore(1000, '\n');

		if (!ifile.good())
			break;
	}

	auto compareFunc = [](const particle &left, const particle &right) -> bool {
		return left.id < right.id;
	};
	std::sort(ps.begin(), ps.end(), compareFunc);

	for (auto p : ps)
		pak.Insert(AtomInfo("A", p.radius), pak.CartesianCoord2RelativeCoord(p.car));
	return pak;
}
Configuration readOmidPacking(const std::string &prefix)
{
	std::fstream ifile(prefix + std::string(".txt"), std::fstream::in);
	return readOmidPacking(ifile);
}

bool file_exist(std::stringstream &fileName)
{
	std::string fn(fileName.str());
	fn += ".txt";
	std::ifstream infile(fn);
	return infile.good();
}

void shear2D(Configuration &c, double strain)
{
	GeometryVector bas[2] = {c.GetBasisVector(0), c.GetBasisVector(1)};
	bas[1].x[0] += strain * bas[0].x[0];
	c.ChangeBasisVector(bas);
}
int Debug()
{
	{
		typedef std::vector<std::pair<int, std::string>> sortable_output;
		double neighbors_rangeCoeff;
		//had to hard-code it because cin cannot tolerate spaces
		//std::string prefix = "/media/scratch/Synced/DEM_data/AR=2.56/ForBetweenness4-long/DEM/postxyz/particle", outFilePrefix;
		std::string prefix = "/mnt/DEMDAA/AR=2.56/ForBetweenness4-long/DEM/postxyz/particle", d2minout = "/mnt/DEMDAA/AR=2.56/ForBetweenness4-long/d2min", outFilePrefix;
		int numType;
		int minTraj, maxTraj, trajIncrement, maxDataSetSize, minusStep, d2min_fileout;
		std::cin >> minTraj >> maxTraj >> trajIncrement >> maxDataSetSize;
		std::cin >> neighbors_rangeCoeff;
		std::cin >> outFilePrefix;
		std::cin >> numType;
		std::vector<std::pair<int, std::stringstream>> vsoft, hard;
		double d2minBound, d2minRange;
		std::cin >> d2minBound >> d2minRange;
		std::cin >> minusStep;
		std::cin >> d2min_fileout;
		std::cerr << minusStep << "\n";

		const double rmin = 0.0045, rmax = 0.0075;
		const double logRRatioMax = std::log(rmax / rmin);

		//need an array of fstreams, which cannot be copied
		//so use unique_ptr instead
		std::vector<std::unique_ptr<std::fstream>> ofiles, ofileh;
		std::vector<sortable_output> ovsoft, ovhard;
		for (int i = 0; i < numType; i++)
		{
			std::stringstream sss, ssh;
			sss << outFilePrefix << "soft" << i << ".txt";
			ssh << outFilePrefix << "hard" << i << ".txt";
			ofiles.push_back(std::move(std::unique_ptr<std::fstream>(new std::fstream(sss.str(), std::fstream::out))));
			ofileh.push_back(std::move(std::unique_ptr<std::fstream>(new std::fstream(ssh.str(), std::fstream::out))));
			ovsoft.push_back(sortable_output(maxDataSetSize));
			ovhard.push_back(sortable_output(maxDataSetSize));
		}

		auto addSfFunc = [&](const Configuration &c0, size_t i, std::ostream &output, const AtomInfo &info) -> void {
			double radius = info.radius;

			std::stringstream ss;
			ss << radius << " ";

			double neighbor_range = 2 * rmax * neighbors_rangeCoeff;

			auto outputFunc = [&](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t Sourceparticle) -> void {
				double r = std::sqrt(shift.Modulus2());

				//new trial : normalize by effective contact radius
				double d1 = radius;
				double d2 = c0.GetCharacteristics(Sourceparticle).radius;
				double d = (d1 + d2);
				double normalizedR = r / d;
				if (r > 0.0 && normalizedR < neighbors_rangeCoeff)
					ss << (char)('A' + std::floor((logRRatioMax + std::log(c0.GetCharacteristics(Sourceparticle).radius / radius)) / (2 * logRRatioMax) * numType)) << " " << normalizedR << " ";
			};
			c0.IterateThroughNeighbors(i, neighbor_range, outputFunc);
			ss << "\n";
#pragma omp critical(output)
			{
				output << ss.str();
				output.flush();
			}
		};

		//New Outputfunction with sortable vectors.
		auto addSfVecFunc = [&](const Configuration &c0, size_t i, sortable_output &output, const AtomInfo &info, const int ts) -> void {
			double radius = info.radius;
			std::stringstream ss;
			ss << ts << " " << radius << " ";

			double neighbor_range = 2 * rmax * neighbors_rangeCoeff;

			auto outputFunc = [&](const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t Sourceparticle) -> void {
				double r = std::sqrt(shift.Modulus2());

				//new trial : normalize by effective contact radius
				double d1 = radius;
				double d2 = c0.GetCharacteristics(Sourceparticle).radius;
				double d = (d1 + d2);
				double normalizedR = r / d;
				if (r > 0.0 && normalizedR < neighbors_rangeCoeff)
					ss << (char)('A' + std::floor((logRRatioMax + std::log(c0.GetCharacteristics(Sourceparticle).radius / radius)) / (2 * logRRatioMax) * numType)) << " " << normalizedR << " ";
			};
			c0.IterateThroughNeighbors(i, neighbor_range, outputFunc);
			ss << "\n";
#pragma omp critical(output)
			{
				output.push_back(std::make_pair(ts, ss.str()));
			}
		};
		std::map<int, std::pair<double, double>> ts_d2min, ts_affine;
		std::vector<size_t> nsoft(numType, 0), nhard(numType, 0);
#pragma omp parallel for
		for (int j = minTraj; j < maxTraj - 1; j += trajIncrement)
		{
			std::stringstream ss0, ss1, ss2;
			ss0 << prefix << j;
			ss1 << prefix << (j + trajIncrement);
			int offset = j - (minusStep * trajIncrement);
			ss2 << prefix << offset;

			/*std::cerr << "ss0: " << ss0.str() << "\n";
			std::cerr  << "ss1: " << ss1.str() << "\n";
			std::cerr  << "ss2: " << ss2.str() << "\n";
	        if (!(file_exist(ss0) && file_exist(ss1) && file_exist(ss2)))
			{
				std::cerr << file_exist(ss0) << " " << ss2.str() << " " << file_exist(ss1) << " " << file_exist(ss2) << "\n";
				continue;
			}*/
			Configuration c0 = readOmidPacking(ss0.str());
			Configuration c1 = readOmidPacking(ss1.str());
			Configuration c2 = readOmidPacking(ss2.str());

			if (c0.NumParticle() != 0)
			{
				std::vector<double> d2min(c0.NumParticle(), 0.0);

				//long double avg = std::accumulate(d2min.begin(), d2min.end(), 0.0);
				if (d2min_fileout)
				{
					std::vector<double> affine_trans(c0.NumParticle(), 0.0);
					for (int i = 0; i < c0.NumParticle(); i++)
						d2min[i] = ::D2Min(c0, c1, i, d2minRange, &affine_trans[i]);
					std::stringstream fname_d2min, fname_affine;
					fname_d2min << d2minout << "/"
								<< "d2min_" << std::setw(9) << std::setfill('0') << j << ".txt";
					fname_affine << d2minout << "/"
								 << "affine_" << std::setw(9) << std::setfill('0') << j << ".txt";
					std::ofstream outf_d2min(fname_d2min.str());
					std::ofstream outf_affine(fname_affine.str());
					double d2min_sum{0};
					for (auto const &e : d2min)
					{
						if (!std::isnan(e))
							d2min_sum += e;
						outf_d2min << e << "\n";
					}
					outf_d2min.close();
					double affine_sum{0};
					for (auto const &e : affine_trans)
					{
						if (!std::isnan(e))
							affine_sum += e;
						outf_affine << e << "\n";
					}
					outf_affine.close();
#pragma omp critical
					{
						ts_d2min[j] = std::make_pair(d2min_sum, d2min_sum / c0.NumParticle());
						ts_affine[j] = std::make_pair(affine_sum, affine_sum / c0.NumParticle());
					}
				}
				else
				{
					for (int i = 0; i < c0.NumParticle(); i++)
						d2min[i] = ::D2Min(c0, c1, i, d2minRange);
				}
				for (int i = 0; i < c2.NumParticle(); i++)
				{
					int type = std::floor(std::log(c2.GetCharacteristics(i).radius / rmin) / logRRatioMax * numType);
					//std::cerr << type << "\n";
					if (type < 0)
					{
						std::cerr << "warning : radius less than rmin, file=" << ss0.str() << ", i=" << i << ", radius=" << c2.GetCharacteristics(i).radius << std::endl;
						continue;
					}
					else if (type >= numType)
					{
						std::cerr << "warning : radius less than rmin, file=" << ss0.str() << ", i=" << i << ", radius=" << c2.GetCharacteristics(i).radius << std::endl;
						continue;
					}
					bool outputS = false, outputH = false;
#pragma omp critical(determineType)
					{
						if (d2min[i] > d2minBound && nsoft[type] < maxDataSetSize)
						{
							nsoft[type]++;
							outputS = true;
						}
						else if (d2min[i] < d2minBound && nhard[type] < nsoft[type])
						{
							nhard[type]++;
							outputH = true;
						}
					}

					//debug temp
					//std::cout<<"radius="<<c0.GetCharacteristics(i).radius<<", type="<<type<<", outputS="<<outputS<<", outputH="<<outputH<<std::endl;
					/*
					if (outputS)
						addSfFunc(c2, i, *ofiles[type], c2.GetCharacteristics(i));
					else if (outputH)
						addSfFunc(c2, i, *ofileh[type], c2.GetCharacteristics(i));
					*/
					if (outputS)
						addSfVecFunc(c2, i, ovsoft[type], c2.GetCharacteristics(i), j);
					else if (outputH)
						addSfVecFunc(c2, i, ovhard[type], c2.GetCharacteristics(i), j);
				}
			}
		}
#pragma omp parallel for
		for (int i = 0; i < numType; i++)
		{
			std::sort(ovhard[i].begin(), ovhard[i].end());
			std::sort(ovsoft[i].begin(), ovsoft[i].end());
			for (const auto &e : ovhard[i])
			{
				*ofileh[i] << e.second;
			}
			for (const auto &e : ovsoft[i])
			{
				*ofiles[i] << e.second;
			}
		}

		if (d2min_fileout)
		{
			std::stringstream ss_d2min, ss_affine;
			ss_d2min << d2minout << "/ts_d2min.txt";
			ss_affine << d2minout << "/ts_affine.txt";
			std::ofstream of_d2min(ss_d2min.str()), of_affine(ss_affine.str());
			of_d2min << "ts sum avg\n";
			of_affine << "ts sum avg\n";
#pragma omp parallel sections
			{
#pragma omp section
				{
					for (auto const &elem : ts_d2min)
					{
						of_d2min << elem.first << " " << elem.second.first << " " << elem.second.second << "\n";
					}
				}
#pragma omp section 
				{
					for (auto const &elem : ts_affine)
					{
						of_affine << elem.first << " " << elem.second.first << " " << elem.second.second << "\n";
					}
				}
			}
		}
		return 0;
	}
}
