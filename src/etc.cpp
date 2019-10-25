#include "etc.h"
#include <cmath>

std::fstream logfile("DetailedLog.txt", std::fstream::app|std::fstream::out);

//unsigned long MCThreadsPerConfiguration=0;
unsigned long MCNumParallelConfigurations=0;
unsigned long PTParallelConfigurationsPerSystem=8;
unsigned long OptimizationNumThreadsPerTrial=0;
unsigned long OptimizationNumThreads=7;

bool AlsoWriteEPS=false;

double TwoBodyDistance_MaxLength=4;

double PlotPhononFrequencies_MaxFrequencySquared_Override=0.0;
size_t Verbosity = 2;
time_t ProgramStart;
time_t TimeLimit;

#include <gsl/gsl_sf_gamma.h>

double HyperSphere_C(DimensionType n)
{
	return std::pow( ::pi, n/2.0)/gsl_sf_gamma(n/2.0+1);
}
double HyperSphere_SurfaceArea(DimensionType n, double R)
{
	return n*HyperSphere_C(n)*std::pow(R, static_cast<double>(n)-1.0);
}
double HyperSphere_Volume(DimensionType n, double R)
{
	return HyperSphere_C(n)*std::pow(R, static_cast<double>(n));
}


#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
//do a linear fit y=c1*x1+c2*x2+...+cn*xn
//input y[i] is y of the ith data point
//      x[j][i] is xj of the ith data point
//return a vector of <c1, c2, ..., cn>
std::vector<double> MultiVariableLinearFit(const std::vector<double> & y, const std::vector< std::vector<double> > & x, double * pchisq, double * pR2)
{
	size_t ndata = y.size();
	size_t nvar = x.size();
	gsl_multifit_linear_workspace * pwor = gsl_multifit_linear_alloc(ndata, nvar);
	gsl_matrix * px = gsl_matrix_alloc(ndata, nvar);
	gsl_matrix * pcov = gsl_matrix_alloc(nvar, nvar);
	gsl_vector * py = gsl_vector_alloc(ndata);
	gsl_vector * pc = gsl_vector_alloc(nvar);

	for (size_t i = 0; i < ndata; i++)
		*gsl_vector_ptr(py, i) = y[i];
	for (size_t i = 0; i < ndata; i++)
		for (size_t j = 0; j < nvar; j++)
			*gsl_matrix_ptr(px, i, j) = x[j][i];

	double chisq;
	gsl_multifit_linear(px, py, pc, pcov, &chisq, pwor);

	std::vector<double> result(nvar, 0.0);
	for (size_t i = 0; i < nvar; i++)
		result[i] = gsl_vector_get(pc, i);

	if (pchisq != nullptr)
		(*pchisq) = chisq;
	if (pR2 != nullptr)
		(*pR2) = 1.0 - chisq / gsl_stats_tss(py->data, py->stride, py->size);

	gsl_vector_free(pc);
	gsl_vector_free(py);
	gsl_matrix_free(pcov);
	gsl_matrix_free(px);
	gsl_multifit_linear_free(pwor);

	return result;
}

#include <chrono>
long GetPreciseClock(void)
{
	std::chrono::system_clock sc;
	return sc.now().time_since_epoch().count();
}
