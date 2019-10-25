#include "StructureFactor.h"

#include <complex>
#include <cmath>
#include <algorithm>

double StructureFactor(const Configuration & Config, const GeometryVector & k)
{
	std::complex<double> rho;
	size_t NumParticle = Config.NumParticle();
	for(size_t i=0; i<NumParticle; i++)
	{
		rho+=std::exp(std::complex<double>(0, k.Dot(Config.GetCartesianCoordinates(i))));
	}
	return (rho.real()*rho.real()+rho.imag()*rho.imag())/NumParticle;
}

std::vector<GeometryVector> GetKs(const PeriodicCellList<Empty> & tempList, double CircularKMax, double LinearKMax, double SampleProbability)
{
	RandomGenerator gen(98765);
	std::vector<GeometryVector> results;
	DimensionType dim = tempList.GetDimension();
	if (dim > 1)
	{
		std::vector<GeometryVector> bs;
		for (DimensionType i = 0; i<dim; i++)
			bs.push_back(tempList.GetReciprocalBasisVector(i));
		//a periodic list of reciprocal lattice
		PeriodicCellList<Empty> reciprocal(dim, &bs[0], std::sqrt(bs[0].Modulus2()));
		reciprocal.Insert(Empty(), GeometryVector(dim));

		std::vector<GeometryVector> ks;
		std::vector<GeometryVector> preKs;
		reciprocal.IterateThroughNeighbors(GeometryVector(dim), CircularKMax, [&ks, &dim](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void{
			bool Trivial = true;
			for (DimensionType d = 0; d<dim; d++)
			{
				if (PeriodicShift[d]>0)
				{
					Trivial = false;
					break;
				}
				else if (PeriodicShift[d]<0)
				{
					Trivial = true;
					break;
				}
			}
			if (!Trivial)
				ks.push_back(shift);
		});
		std::sort(ks.begin(), ks.end(), [](const GeometryVector & left, const GeometryVector & right) ->bool {return left.Modulus2()<right.Modulus2(); });
		//ks contains K points in that circle

		if (LinearKMax == CircularKMax)
			return ks;

		std::vector<GeometryVector> tempBase;
		for (DimensionType i = 0; i<dim; i++)
		{
			GeometryVector t(dim);
			t.x[i] = 2;
			tempBase.push_back(t);
		}
		PeriodicCellList<Empty> KDirection(dim, &tempBase[0], std::sqrt(tempBase[0].Modulus2())*std::pow(1.0 / ks.size(), 1.0 / (dim))*10.0);//use this list to discover k points of the same direction
		for (auto iter = ks.begin(); iter != ks.end(); iter++)
		{
			double Length = std::sqrt(iter->Modulus2());
			if (Length == 0)
				continue;
			GeometryVector dir(dim);
			for (DimensionType i = 0; i<dim; i++)
				dir.x[i] = iter->x[i] / Length / 2.0;

			bool NeighborFound = false;
			KDirection.IterateThroughNeighbors(dir, 1e-10, [&NeighborFound](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void{NeighborFound = true; });
			if (NeighborFound == false)
			{
				preKs.push_back(*iter);
				KDirection.Insert(Empty(), dir);
			}
		}
		//preKs contains K points in that circle with different directions

		ks.clear();
		for (auto iter = preKs.begin(); iter != preKs.end(); iter++)
		{
			for (size_t i = 1;; i++)
			{
				GeometryVector temp = static_cast<double>(i)*(*iter);
				if (temp.Modulus2()>LinearKMax*LinearKMax)
					break;
				else
					if(SampleProbability>=1.0 || gen.RandomDouble()<SampleProbability)
						results.push_back(temp);
			}
		}
		preKs.clear();
	}
	else
	{
		GeometryVector b = tempList.GetReciprocalBasisVector(0);
		for (double i = 1.0;; i += 1.0)
		{
			GeometryVector k = i*b;
			if (k.x[0] > LinearKMax)
				break;
			if (SampleProbability >= 1.0 || gen.RandomDouble()<SampleProbability)
				results.push_back(k);
		}
	}
	return results;
}


//calculate some function of some arbitrary type ConfigType and the reciprocal lattice vector k, F(k), averag over directions and all results from GetConfigsFunction, and bin them
//Fk is a function that takes a ConfigType and a k vector and produces F(k)
//ConfigType must support all operations of a PeriodicCellList
template<typename ConfigType, typename FkType> void IsotropicFk(FkType Fk, std::function<const ConfigType(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision, double SampleProbability)
{
	Results.clear();
	if(!(CircularKMax>0))
		return;
	if(!(LinearKMax>0))
		return;

	struct SkBin
	{
		double Sum1, SumK2, SumS, SumS2;
		SkBin() : Sum1(0), SumK2(0), SumS(0), SumS2(0)
		{}
	};
	if (KPrecision == 0.0)
	{
		std::cerr << "Warning in IsotropicStructureFactor : This version does not support KPrecision==0.0. Auto choosing this quantity!\n";
		ConfigType c = GetConfigsFunction(0);
		KPrecision = std::sqrt(c.GetReciprocalBasisVector(0).Modulus2());
	}
	size_t NumBin = std::floor(LinearKMax / KPrecision) + 1;
	std::vector<SkBin> vSkBin(NumBin, SkBin());


	GeometryVector prevBasis [ ::MaxDimension];
	std::vector<GeometryVector> ks;
	if(Verbosity>1)
		std::cout<<"Computing S(k)";
	progress_display pd(NumConfigs);
	for(size_t j=0; j<NumConfigs; j++)
	{
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		ConfigType CurrentConfig = GetConfigsFunction(j);
		if(CurrentConfig.GetDimension()==0)
			break;
		if(j!=0)
		{
			bool SameBasis = true;
			for(DimensionType i=0; i<CurrentConfig.GetDimension(); i++)
				if( !(prevBasis[i]==CurrentConfig.GetBasisVector(i)) )
					SameBasis=false;

			if(SameBasis==false)
				ks=GetKs(CurrentConfig, CircularKMax, LinearKMax, SampleProbability);
		}
		else
			ks = GetKs(CurrentConfig, CircularKMax, LinearKMax, SampleProbability);

		signed long end = ks.size();
#pragma omp parallel for schedule(guided)
		for(signed long i=0; i<end; i++)
		{
			double s=Fk(CurrentConfig, ks[i]);
			double k2 = ks[i].Modulus2();
			size_t Bin = std::floor(std::sqrt(k2) / KPrecision);
#pragma omp atomic
			vSkBin[Bin].Sum1 += 1.0;
#pragma omp atomic
			vSkBin[Bin].SumK2 += k2;
#pragma omp atomic
			vSkBin[Bin].SumS += s;
#pragma omp atomic
			vSkBin[Bin].SumS2 += s*s;
		}
		for(DimensionType i=0; i<CurrentConfig.GetDimension(); i++)
			prevBasis[i]=CurrentConfig.GetBasisVector(i);
		pd++;
	}
	for (auto iter = vSkBin.begin(); iter != vSkBin.end(); iter++)
	{
		if (iter->Sum1 != 0.0)
		{
			GeometryVector temp(4);
			temp.x[0] = std::sqrt(iter->SumK2 / iter->Sum1);
			temp.x[1] = iter->SumS / iter->Sum1;
			temp.x[2] = KPrecision;
			temp.x[3] = std::sqrt((iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1));
			Results.push_back(temp);
		}
	}
	
	if(Verbosity>2)
		std::cout<<"done!\n";
}


void IsotropicStructureFactor(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision, double SampleProbability)
{
	::IsotropicFk(StructureFactor, GetConfigsFunction, NumConfigs, CircularKMax, LinearKMax, Results, KPrecision, SampleProbability);
}


#include <boost/math/special_functions.hpp>
double Mtilde(double k, double R, DimensionType d)
{
	if (d == 1)
		return 2 * std::sin(k*R) / k;
	else if (d == 2)
		return 2 * pi*R*boost::math::cyl_bessel_j(1, k*R) / k;
	else if (d == 3)
		return 4 * pi*(std::sin(k*R) - k*R*std::cos(k*R)) / k / k / k;
	else if (d == 4)
		return 4 * pi*pi*R*(2 * boost::math::cyl_bessel_j(1, k*R) - k*R*boost::math::cyl_bessel_j(0, k*R)) / k / k / k;
	else
	{
		std::cerr << "Error in Mtilde : unsupported dimension!\n";
		return 0.0;
	}
}

double SpectralDensity(const SpherePacking & Config, const GeometryVector & k)
{
	std::complex<double> rho;
	std::map<double, double> mtildes;
	size_t NumParticle = Config.NumParticle();
	double kl = std::sqrt(k.Modulus2());
	DimensionType d = Config.GetDimension();
	for (size_t i = 0; i<NumParticle; i++)
	{
		std::complex<double> a = std::exp(std::complex<double>(0, k.Dot(Config.GetCartesianCoordinates(i))));
		double mt;
		double R = Config.GetCharacteristics(i);
		auto iter = mtildes.find(R);
		if (iter != mtildes.end())
			mt = iter->second;
		else
		{
			mt = Mtilde(kl, R, d);
			mtildes.insert(std::make_pair(R, mt));
		}
		rho += a*mt;
	}
	//debug temp
	//std::cout << rho << '\n';
	return (rho.real()*rho.real() + rho.imag()*rho.imag()) / Config.PeriodicVolume();
}

//todo : test this function
void IsotropicSpectralDensity(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision, double SampleProbability)
{
	::IsotropicFk(SpectralDensity, GetConfigsFunction, NumConfigs, CircularKMax, LinearKMax, Results, KPrecision, SampleProbability); 
}

//old code

//void IsotropicSpectralDensity(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision)
//{
//	Results.clear();
//	if (!(CircularKMax>0))
//		return;
//	if (!(LinearKMax>0))
//		return;
//
//	struct SkBin
//	{
//		double Sum1, SumK2, SumS, SumS2;
//		SkBin() : Sum1(0), SumK2(0), SumS(0), SumS2(0)
//		{}
//	};
//	if (KPrecision == 0.0)
//	{
//		std::cerr << "Warning in IsotropicSpectralDensity : This version does not support KPrecision==0.0. Auto choosing this quantity!\n";
//		SpherePacking c = GetConfigsFunction(0);
//		KPrecision = std::sqrt(c.GetReciprocalBasisVector(0).Modulus2());
//	}
//	size_t NumBin = std::floor(LinearKMax / KPrecision) + 1;
//	std::vector<SkBin> vSkBin(NumBin, SkBin());
//
//
//	GeometryVector prevBasis[::MaxDimension];
//	std::vector<GeometryVector> ks;
//	if (Verbosity>1)
//		std::cout << "Computing chi(k)";
//	progress_display pd(NumConfigs);
//	for (size_t j = 0; j<NumConfigs; j++)
//	{
//		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
//		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
//		SpherePacking CurrentConfig = GetConfigsFunction(j);
//		if (CurrentConfig.GetDimension() == 0)
//			continue;
//		if (j != 0)
//		{
//			bool SameBasis = true;
//			for (DimensionType i = 0; i<CurrentConfig.GetDimension(); i++)
//				if (!(prevBasis[i] == CurrentConfig.GetBasisVector(i)))
//					SameBasis = false;
//
//			if (SameBasis == false)
//				ks=GetKs(CurrentConfig, CircularKMax, LinearKMax);
//		}
//		else
//			ks=GetKs(CurrentConfig, CircularKMax, LinearKMax);
//
//		signed long end = ks.size();
//#pragma omp parallel for schedule(guided)
//		for (signed long i = 0; i<end; i++)
//		{
//			double s = SpectralDensity(CurrentConfig, ks[i]);
//			double k2 = ks[i].Modulus2();
//			size_t Bin = std::floor(std::sqrt(k2) / KPrecision);
//#pragma omp atomic
//			vSkBin[Bin].Sum1 += 1.0;
//#pragma omp atomic
//			vSkBin[Bin].SumK2 += k2;
//#pragma omp atomic
//			vSkBin[Bin].SumS += s;
//#pragma omp atomic
//			vSkBin[Bin].SumS2 += s*s;
//		}
//		for (DimensionType i = 0; i<CurrentConfig.GetDimension(); i++)
//			prevBasis[i] = CurrentConfig.GetBasisVector(i);
//		pd++;
//	}
//	for (auto iter = vSkBin.begin(); iter != vSkBin.end(); iter++)
//	{
//		if (iter->Sum1 != 0.0)
//		{
//			GeometryVector temp(4);
//			temp.x[0] = std::sqrt(iter->SumK2 / iter->Sum1);
//			temp.x[1] = iter->SumS / iter->Sum1;
//			temp.x[2] = KPrecision;
//			temp.x[3] = std::sqrt(iter->SumS2 / (iter->Sum1) - temp.x[1] * temp.x[1]) / (iter->Sum1);
//			Results.push_back(temp);
//		}
//	}
//
//	if (Verbosity>2)
//		std::cout << "done!\n";
//}
