#ifndef ANALYZEDATA_INCLUDED
#define ANALYZEDATA_INCLUDED

#include <vector>

const double SampleAutocorrelationLimit=0.01;//when calculating the error, we need to integrate autocorrelation, when abs(autocorrelation) is smaller than this value, cutoff the integral

inline void AnalyzeData(const std::vector<double> & data, double & average, double & error, double & corrtime)
{
	//calculate <data> and <data^2>
	double a2=0;
	average=0;
	for(auto iter=data.begin(); iter!=data.end(); iter++)
	{
		average+=*iter;
		a2+=(*iter) * (*iter);
	}
	average/=data.size();
	a2/=data.size();

	//calculate autocorrelation function
	corrtime=0.5;
	for(int i=0; i<data.size(); i++)
	{
		double autocor=0.0;
		for(auto iter=data.begin()+i; iter!=data.end(); iter++)
			autocor+=(*(iter-i)-average)*(*(iter)-average);
		autocor/=(data.size()-i);
		autocor=(autocor)/(a2-average*average);

		corrtime+=autocor;

		if(std::abs(autocor)<SampleAutocorrelationLimit) break;
		if(autocor!=autocor) break;//this is true if and only if autocor=NaN
	}

	error=std::sqrt((a2-average*average)*(1+2*corrtime)/data.size());
}

inline void AnalyzeData(const std::vector<double> & data, double & mean, double & stddev)
{
	double sum = 0.0, sum2 = 0.0;
	for (auto iter = data.begin(); iter != data.end(); iter++)
	{
		sum += *iter;
		sum2 += (*iter)*(*iter);
	}
	mean = sum / data.size();
	stddev = std::sqrt(sum2 / data.size() - mean*mean);
}


#endif