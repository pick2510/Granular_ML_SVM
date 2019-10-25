#ifndef PLOTS_INCLUDED
#define PLOTS_INCLUDED

#include "ParameterSet.h"
#include "LatticeSumSolver.h"


//parameter ParticleSizeRescale only works for 2D configurations
void Plot(std::string prefix, const Configuration & c, const std::string & remark="", double ParticleSizeRescale=1.0);

void Plots(IsotropicPairPotential * pPot, const std::vector<LatticeSumSolver *> * Solvers, const char * Prefix, size_t Enphasis=0);
void Plots(const ParameterSet & Param, const std::vector<LatticeSumSolver *> * Solvers, const char * Prefix, size_t Enphasis=0);

//require that all configurations in this vector have the same basis vectors
void PlotDirectionalStructureFactor2D(const std::vector<Configuration> & vc, long limit, std::string prefix, bool LogScale, double MinSOverride, double MaxSOverride, double TileSize);
void PlotDirectionalStructureFactor2D(const Configuration & c, long limit, std::string prefix, bool LogScale=true, double MinSOverride=0.0, double MaxSOverride=0.0, double TileSize=1.0);
inline void PlotDirectionalStructureFactor2D(std::string prefix, const Configuration & c, long limit, bool LogScale=true, double MinSOverride=0.0, double MaxSOverride=0.0, double TileSize=1.0)
{
	PlotDirectionalStructureFactor2D(c, limit, prefix, LogScale, MinSOverride, MaxSOverride, TileSize);
}

//Plot a path in a 2D configuration
void PlotPath(std::string prefix, const Configuration & List, const std::string & remark, const std::vector<size_t> & path);

void PlotDisplacement(std::string prefix, const Configuration & List, const std::string & remark, const std::vector<double> & displacement);
template<typename add>
void PlotPacking2D(std::string prefix, const PeriodicCellList<add> & List, const std::string & remark, double r);
void PlotPacking2D(std::string prefix, const SpherePacking & List, const std::string & remark = "", std::function<std::string(size_t Num)> GetMGLOptionsFunc = [](size_t Num){return "b"; });


void PlotFunction_MathGL(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel);


//functions and classes to output xmGrace/qtGrace figures
//Use the class if you want even more control
extern std::string PlotFunction_Grace_SetType;
void PlotFunction_Grace(const std::vector<GeometryVector> * presult, size_t NumSet, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel, const std::vector<std::string> & legends, const std::string & Title, double MinX, double MaxX, double TickX, double MinY, double MaxY, double TickY);
void PlotFunction_Grace(const std::vector<GeometryVector> * presult, size_t NumSet, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel, const std::vector<std::string> & legends, const std::string & Title);
void PlotFunction_Grace(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix="temp", const std::string & xLabel="", const std::string & yLabel="", const std::string & Title="");


class gracePlot
{
	const int NColor = 14;
	const int NLineStyle = 8;
	const int NSymbol = 10;
public:
	class dataSet
	{
	public:
		std::vector<GeometryVector> result;
		std::string type;//xy, xydy, or xydxdy
		int lineColor;
		int lineStyle;
		int symbolColor;
		int symbolStyle;
		std::string legend;
		double symbolSize;
	};
	std::vector<dataSet> dataSets;
	double MinX, MaxX, TickX, MinY, MaxY, TickY;
	std::string xLabel, yLabel, Title;
	std::string xScale, yScale;//"Normal" or "Logarithmic"
	gracePlot() : MinX(0.0), MaxX(1.0), TickX(0.5), MinY(0.0), MaxY(1.0), TickY(0.5), xScale("Normal"), yScale("Normal")
	{}
	void autoScale();
	void autoTick();

	void autoScaleAndTick();
	void addDataSet(std::vector<GeometryVector> data, std::string legend = "", int lineColor = -1, int lineStyle = -1, int symbolColor = -1, int symbolStyle = -1, std::string type = "xy", double symbolSize = 1.0); // for colors and styles, -1 means auto, and 0 means none
	void outputFigure(const std::string & OutFilePrefix);
};


//ToDo : revamp them as constructors of class GracePlot
void ReadGraceData(std::vector<std::vector<GeometryVector> > & result, const std::string & InFilePrefix);
void ReadGraceData(std::vector<std::vector<GeometryVector> > & result, std::istream & ifile);

void PlotG1_2D(const std::vector<Configuration> & vc, long numBinsPerSide, std::string prefix, std::string remark);
#endif