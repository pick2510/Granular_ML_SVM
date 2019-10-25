#ifndef ETC_INCLUDED
#define ETC_INCLUDED


const double LengthPrecision=1e-4;
const double EnergyPrecision=1e-4;
const double MaxDistance = 1e10;
const double MaxEnergy=1e10;
const double pi=3.1415926535897932384626433832795028841971693993751058209;
const double E = 2.71828182845904523536;
typedef unsigned short DimensionType;
const DimensionType MaxDimension = 4;
const bool UseCorrectedPotential = false;

extern bool AlsoWriteEPS;//when write bmp files, also write eps

//extern unsigned long MCThreadsPerConfiguration;
extern unsigned long MCNumParallelConfigurations;
extern unsigned long PTParallelConfigurationsPerSystem;
extern unsigned long OptimizationNumThreadsPerTrial;
extern unsigned long OptimizationNumThreads;

extern double TwoBodyDistance_MaxLength;

extern double PlotPhononFrequencies_MaxFrequencySquared_Override;

#ifdef USE_PNG
const char * const FigureFormat = ".png";
#else
const char * const FigureFormat = ".bmp";
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

extern std::fstream logfile;

extern size_t Verbosity;

struct Empty
{
	Empty()
	{}
};

extern time_t ProgramStart;

//implemented in : structure optimization, 
extern time_t TimeLimit;


double HyperSphere_C(DimensionType n);
double HyperSphere_SurfaceArea(DimensionType n, double R);
double HyperSphere_Volume(DimensionType n, double R);


class progress_display
{
public:
	explicit progress_display( unsigned long expected_count,
		std::ostream & os = (Verbosity>1) ? std::cout : logfile ,
		const std::string & s1 = "\n", //leading strings
		const std::string & s2 = "",
		const std::string & s3 = "" )
		// os is hint; implementation may ignore, particularly in embedded systems
		: m_os(os), m_s1(s1), m_s2(s2), m_s3(s3) { restart(expected_count); }

	void restart( unsigned long expected_count )
		//  Effects: display appropriate scale
		//  Postconditions: count()==0, expected_count()==expected_count
	{
		_count = _next_tic_count = _tic = 0;
		_expected_count = expected_count;

		m_os << m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
			<< m_s2 << "|----|----|----|----|----|----|----|----|----|----|"
			<< std::endl  // endl implies flush, which ensures display
			<< m_s3;
		if ( !_expected_count ) _expected_count = 1;  // prevent divide by zero
	} // restart

	unsigned long  operator+=( unsigned long increment )
		//  Effects: Display appropriate progress tic if needed.
		//  Postconditions: count()== original count() + increment
		//  Returns: count().
	{
		if ( (_count += increment) >= _next_tic_count ) { display_tic(); }
		return _count;
	}

	unsigned long  operator++()           { return operator+=( 1 ); }
	unsigned long  operator++(int)           { return (operator+=( 1 )-1); }
	unsigned long  count() const          { return _count; }
	unsigned long  expected_count() const { return _expected_count; }

private:
	std::ostream &     m_os;  // may not be present in all imps
	const std::string  m_s1;  // string is more general, safer than 
	const std::string  m_s2;  //  const char *, and efficiency or size are
	const std::string  m_s3;  //  not issues

	unsigned long _count, _expected_count, _next_tic_count;
	unsigned int  _tic;
	void display_tic()
	{
		// use of floating point ensures that both large and small counts
		// work correctly.  static_cast<>() is also used several places
		// to suppress spurious compiler warnings. 
		unsigned int tics_needed =
			static_cast<unsigned int>(
			(static_cast<double>(_count)/_expected_count)*50.0 );
		do { m_os << '*' << std::flush; } while ( ++_tic < tics_needed );
		_next_tic_count = 
			static_cast<unsigned long>((_tic/50.0)*_expected_count);
		if ( _count == _expected_count ) {
			if ( _tic < 51 ) m_os << '*';
			m_os << std::endl;
		}
	} // display_tic
};

#include <vector>
//do a linear fit y=c1*x1+c2*x2+...+cn*xn
//input y[i] is y of the ith data point
//      x[j][i] is xj of the ith data point
//return a vector of <c1, c2, ..., cn>
std::vector<double> MultiVariableLinearFit(const std::vector<double> & y, const std::vector< std::vector<double> > & x, double * pchisq = nullptr, double * pR2 = nullptr);

//In windows, the return value of this function changes 10,000,000 every second
long GetPreciseClock(void);

struct KahanAccumulation
{
	double sum;
	double correction;
	KahanAccumulation() : sum(0.0), correction(0.0) {}
};

inline KahanAccumulation KahanSum(KahanAccumulation accumulation, double value)
{
	KahanAccumulation result;
	double y = value - accumulation.correction;
	double t = accumulation.sum + y;
	result.correction = (t - accumulation.sum) - y;
	result.sum = t;
	return result;
}


#endif
