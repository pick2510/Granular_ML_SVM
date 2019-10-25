#ifndef STRUCTUREFACTOR_INCLUDED
#define STRUCTUREFACTOR_INCLUDED

#include "GeometryVector.h"
#include "PeriodicCellList.h"
#include <functional>
#include <vector>


//find all ks in the range CircularKMax (i.e. for all ks that abs(k)<CircularKMax) and all of their multiples in the range LinearKMax
std::vector<GeometryVector> GetKs(const PeriodicCellList<Empty> & tempList, double CircularKMax, double LinearKMax, double SampleProbability = 1.0);

//calculate the Structure factor of certain k
double StructureFactor(const Configuration & Config, const GeometryVector & k);

//calculate structure factor of all ks in the range CircularKMax (i.e. for all ks that abs(k)<CircularKMax) and all of their multiples in the range LinearKMax
//Results has the following form:
//( abs(k), S(k), KPrecision, \delta S(k) )
//GetConfigsFunction should return configurations of index i when called
//if SampleProbability is <1, then for any k point satisfying the above condition there is SampleProbability probability that it will be used to calculate S(k) and (1-SampleProbability) probability that it will not be used.
void IsotropicStructureFactor(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision=0.01, double SampleProbability=1.0);


//fourier transform of the indicator function of a hypersphere of radius R
double Mtilde(double k, double R, DimensionType d);


double SpectralDensity(const SpherePacking & Config, const GeometryVector & k);

//todo : test this function
void IsotropicSpectralDensity(std::function<const SpherePacking(size_t i)> GetConfigsFunction, size_t NumConfigs, double CircularKMax, double LinearKMax, std::vector<GeometryVector> & Results, double KPrecision = 0.01, double SampleProbability = 1.0);

//todo : create a class DigitizedConfiguration which supports variable basis vector
//write an overload of IsotropicSpectralDensity for this type, using DigitizedSkCalculator and IsotropicFk
//the first template parameter of IsotropicFk, ConfigType, will need to be a class that includes a DigitizedSkCalculator and the basis vectors
//for the basis vectors part, separate PeriodicCellList into a class SimulationBox and a PeriodicCellList that inherits from it.
#include "KissFFT/kiss_fft.h"
#include "KissFFT/kiss_fftnd.h"
//FFT-based structure factor and spectral density calculator for lattice gases, digitized two-phase media, etc.
//todo : test if this class works well when numPixelPerSide is different for different directions
class DigitizedSkCalculator
{
protected:
	kiss_fft_cpx *inbuf, *outbuf;
	size_t numPixels;
	bool FFTFinished;
	std::vector<int> _numPixelPerSide;
	DimensionType _dim;
public:
	DigitizedSkCalculator(DimensionType dim, const std::vector<int> & numPixelPerSide)
	{
		_dim = dim;
		_numPixelPerSide = numPixelPerSide;
		numPixels = 1;
		for (int i = 0; i < dim; i++)
		{
			if (numPixelPerSide[i] <= 0)
			{
				std::cout << "Error in DigitizedSkCalculator constructor : non-positive numPixelPerSide.\n";
				return;
			}
			numPixels *= numPixelPerSide[i];
		}

		inbuf = new kiss_fft_cpx[numPixels];
		outbuf = new kiss_fft_cpx[numPixels];

		FFTFinished = false;
	}
	// set real-space configuration
	//VoxelType can be bool, char, double, etc., depending on the need.
	//data is a vector of numPixelPerSide[0]*numPixelPerSide[1]*...*numPixelPerSide[dim]
	template <typename PixelType> void SetConfiguration(const std::vector<PixelType> & data)
	{
		for (size_t i = 0; i < numPixels; i++)
		{
			inbuf[i].r = data[i];
			inbuf[i].i = 0.0;
		}
		FFTFinished = false;
	}
	void DoFFT()
	{
		kiss_fftnd_cfg st = kiss_fftnd_alloc(_numPixelPerSide.data(), _dim, false, 0, 0);
		kiss_fftnd(st, inbuf, outbuf);
		free(st);
		FFTFinished = true;
	}
	double GetSk(const std::vector<int> & k, size_t N)
	{
		if (!FFTFinished)
			DoFFT();
		size_t index = 0;
		for (DimensionType d = 0; d < _dim; d++)
		{
			index *= _numPixelPerSide[d];
			index += k[d] % _numPixelPerSide[d];
		}
		kiss_fft_cpx & temp = outbuf[index];
		double s = (temp.i*temp.i + temp.r*temp.r) / N;
		return s;
	}
	double GetChik(const std::vector<int> & k, double Volume)
	{
		return GetSk(k, numPixels*numPixels)*Volume;
	}
	~DigitizedSkCalculator()
	{
		delete[] inbuf;
		delete[] outbuf;
	}
};


#endif