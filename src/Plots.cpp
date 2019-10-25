#include "Plots.h"
#include "Potential.h"
#include <string>
#include <fstream>

#ifdef USE_MATHGL
#include <mgl2/mgl.h>
#endif

#ifdef USE_MATHGL
void PlotCircle2D(mglGraph & graph, double x, double y, double r, const char * MGLOptions, size_t NumLineSegments=100)
{
	mglData xc(NumLineSegments+1), yc(NumLineSegments+1);
	for(size_t i=0; i<NumLineSegments; i++)
	{
		xc.a[i]=x+r*std::cos((double)(i)/NumLineSegments*2* ::pi);
		yc.a[i]=y+r*std::sin((double)(i)/NumLineSegments*2* ::pi);
	}
	xc.a[NumLineSegments]=x+r;
	yc.a[NumLineSegments]=y;
	graph.Plot(xc, yc, MGLOptions);
}
void PlotDisk2D(mglGraph & graph, double x, double y, double r, const char * MGLOptions, size_t NumLineSegments=200)
{
	//size_t HalfNumLineSegments=NumLineSegments/2+1;
	//mglData xc(HalfNumLineSegments+1), yc(HalfNumLineSegments+1), yc2(HalfNumLineSegments+1);
	//for(size_t i=0; i<HalfNumLineSegments; i++)
	//{
	//	xc.a[i]=x+r*std::cos((double)(i)/HalfNumLineSegments* ::pi);
	//	yc.a[i]=y-r*std::sin((double)(i)/HalfNumLineSegments* ::pi);
	//	yc2.a[i]=y+r*std::sin((double)(i)/HalfNumLineSegments* ::pi);
	//}
	//xc.a[HalfNumLineSegments]=x-r;
	//yc.a[HalfNumLineSegments]=y;
	//yc2.a[HalfNumLineSegments]=y;
	//graph.Region(xc, yc, yc2, MGLOptions);
	graph.Circle(mglPoint(x, y, 0.0), r, MGLOptions);
}
#endif


#ifdef USE_MATHGL
void PlotFunction_MathGL(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel)
{
	const size_t NumData=result.size();
	if(NumData==0)
	{
		std::cerr<<"Error in PlotFunction_MathGL: No Data ! Exiting function \n";
		return;
	}
	mglData x(NumData), y(NumData), ex(NumData), ey(NumData);
	double MinX=result[0].x[0]-result[0].x[2], MaxX=result[0].x[0]+result[0].x[2];
	double MinY=result[0].x[1]-result[0].x[3], MaxY=result[0].x[1]+result[0].x[3];
	for(size_t i=0; i<NumData; i++)
	{
		x.a[i]=result[i].x[0];
		y.a[i]=result[i].x[1];
		ex.a[i]=result[i].x[2];
		ey.a[i]=result[i].x[3];
		if(MinX>result[i].x[0]-result[i].x[2])
			MinX=result[i].x[0]-result[i].x[2];
		if(MaxX<result[i].x[0]+result[i].x[2])
			MaxX=result[i].x[0]+result[i].x[2];
		if(MinY>result[i].x[1]-result[i].x[3])
			MinY=result[i].x[1]-result[i].x[3];
		if(MaxY<result[i].x[1]+result[i].x[3])
			MaxY=result[i].x[1]+result[i].x[3];
	}
	mglGraph gr;
	gr.SetRanges(MinX*0.999999999999999, MaxX*1.000000000000001, MinY, MaxY);
	gr.SetOrigin(0, MinY);
	gr.SetFontSize(6);
	gr.SetTickTempl('y', " ");
	//gr.SetTickShift(mglPoint(0.1, 0, 0));
	gr.Axis();
	//gr.SetTickShift(mglPoint(0.0, 0.0));
	gr.SetTickTempl('y', "");
	gr.SetTickTempl('x', " ");
	gr.Axis();
	gr.Label('y', yLabel.c_str(), 0);
	gr.Label('x', xLabel.c_str(), 0);
	gr.SetRanges(MinX, MaxX, MinY, MaxY);
	gr.SetOrigin(MaxX, MaxY);
	gr.SetTickTempl('x', " ");
	gr.SetTickTempl('y', " ");
	gr.Axis();
	//gr.Error(x, y, ex, ey, "k.");
	gr.Plot(x, y, "b .");
	std::string filename;
	filename.assign(OutFilePrefix);
	filename += FigureFormat;
	gr.WriteFrame(filename.c_str());
	if(AlsoWriteEPS)
	{
		filename.assign(OutFilePrefix);
		filename+=".eps";
		gr.WriteEPS(filename.c_str());
	}
}
#else
void PlotFunction_MathGL(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel)
{
	logfile<<"MathGL not Enabled, skip PlotFunction_MathGL!\n";
	return;
}
#endif

void Plot(std::string prefix, const Configuration & List, const std::string & remark, double ParticleSizeRescale)
{
	const Configuration * list= & List;
#ifdef USE_MATHGL

	if(list->GetDimension()==1)
	{
		//calculate the range of the plot
		double max=list->GetBasisVector(0).x[0];
		GeometryVector temp(max, 0.0);
		GeometryVector a0(0, 0);
		GeometryVector a1=a0+temp;
		double height = a1.x[0] / 19.20;
		GeometryVector a2=a0+GeometryVector(0.0, height);
		GeometryVector a3=a2+temp;
		//prepare the graph
		mglGraph gr;
		gr.SetRanges(0.0, max, 0.0, max/19.2);
		gr.Aspect(1.0, 1.0);
		gr.SetSize(1920, 300);
		//gr.SetOrigin(0, 0);
		gr.Puts(mglPoint(0.5*(max), 1.1*(max / 19.2)), remark.c_str(), ":C", 10.0);

		//plot crystal BasisVectors
		{
			mglData x(2), y(2);
			x.a[0]=a0.x[0];
			y.a[0]=a0.x[1];
			x.a[1]=a1.x[0];
			y.a[1]=a1.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a1.x[0];
			y.a[0]=a1.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a0.x[0];
			y.a[1]=a0.x[1];
			gr.Plot(x, y, "k-");
		}
		//plot particles
		double radius=0.1*(max)/list->NumParticle()*ParticleSizeRescale;
		for(int i=0; i<list->NumParticle(); i++)
		{
			GeometryVector now=list->GetCartesianCoordinates(i);
			//PlotDisk2D(gr, now.x[0], 0.5*a2.x[1], radius, "b", 20);
			if(list->GetCharacteristics(i)==AtomInfo("B"))
				gr.Line(mglPoint(now.x[0], 0.3*a2.x[1]), mglPoint(now.x[0], 0.7*a2.x[1]), "q");
			else
				gr.Line(mglPoint(now.x[0], 0.3*a2.x[1]), mglPoint(now.x[0], 0.7*a2.x[1]));
		}

		prefix+=std::string(FigureFormat);
		gr.WriteFrame(prefix.c_str());
		if(AlsoWriteEPS)
		{
			prefix.resize(prefix.size()-4);
			prefix+=std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else if(list->GetDimension()==2)
	{
		//calculate the range of the plot
		GeometryVector a0(0, 0);
		GeometryVector a1=a0+list->GetBasisVector(0);
		GeometryVector a2=a0+list->GetBasisVector(1);
		GeometryVector a3=a2+list->GetBasisVector(0);
		double min=0, max=0;
		if(min>a0.x[0])
			min=a0.x[0];
		if(min>a0.x[1])
			min=a0.x[1];
		if(max<a0.x[0])
			max=a0.x[0];
		if(max<a0.x[1])
			max=a0.x[1];
		if(min>a1.x[0])
			min=a1.x[0];
		if(min>a1.x[1])
			min=a1.x[1];
		if(max<a1.x[0])
			max=a1.x[0];
		if(max<a1.x[1])
			max=a1.x[1];
		if(min>a2.x[0])
			min=a2.x[0];
		if(min>a2.x[1])
			min=a2.x[1];
		if(max<a2.x[0])
			max=a2.x[0];
		if(max<a2.x[1])
			max=a2.x[1];
		if(min>a3.x[0])
			min=a3.x[0];
		if(min>a3.x[1])
			min=a3.x[1];
		if(max<a3.x[0])
			max=a3.x[0];
		if(max<a3.x[1])
			max=a3.x[1];

		//prepare the graph
		mglGraph gr;
		gr.SetRanges(min, max, min, max);
		gr.Aspect(1.0, 1.0);
		gr.SetSize(3200, 3200);
		//gr.SetOrigin(0, 0);
		gr.Puts(mglPoint(0.5*(max+min), min+1.1*(max-min)), remark.c_str());


		//plot crystal BasisVectors
		{
			mglData x(2), y(2);
			x.a[0]=a0.x[0];
			y.a[0]=a0.x[1];
			x.a[1]=a1.x[0];
			y.a[1]=a1.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a1.x[0];
			y.a[0]=a1.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a0.x[0];
			y.a[1]=a0.x[1];
			gr.Plot(x, y, "k-");
		}
		//plot particles
		double radius=0.1*(max-min)/std::sqrt(list->NumParticle())*ParticleSizeRescale;
		//double radius = 0.5;

		for(int i=0; i<list->NumParticle(); i++)
		{
			GeometryVector now=list->GetCartesianCoordinates(i);
			const char * format = "b";
			//debug temp
			if (list->GetCharacteristics(i).name[0] == 'Z')
				format = "r";

			PlotDisk2D(gr, now.x[0], now.x[1], radius, "b", 20);
			//if (list->GetCharacteristics(i) == AtomInfo("A"))
			//	PlotDisk2D(gr, now.x[0], now.x[1], 0.5, "b", 20);
			//else
			//	PlotDisk2D(gr, now.x[0], now.x[1], 0.7, "r", 20);
		}

		prefix += std::string(FigureFormat);
		gr.WriteFrame(prefix.c_str());
		if(AlsoWriteEPS)
		{
			prefix.resize(prefix.size()-4);
			prefix+=std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else if(list->GetDimension()==3)
	{
		//calculate the range of the plot
		std::vector<GeometryVector> a(8, GeometryVector(0.0, 0.0, 0.0));
		a[1]=a[0]+list->GetBasisVector(0);
		a[2]=a[0]+list->GetBasisVector(1);
		a[3]=a[0]+list->GetBasisVector(2);
		a[4]=a[1]+list->GetBasisVector(1);
		a[5]=a[1]+list->GetBasisVector(2);
		a[6]=a[2]+list->GetBasisVector(2);
		a[7]=a[4]+list->GetBasisVector(2);
		double min=0, max=0;
		for(int i=0; i<8; i++)
		{
			if(min>a[i].x[0])
				min=a[i].x[0];
			if(min>a[i].x[1])
				min=a[i].x[1];
			if(min>a[i].x[2])
				min=a[i].x[2];
			if(max<a[i].x[0])
				max=a[i].x[0];
			if(max<a[i].x[1])
				max=a[i].x[1];
			if(max<a[i].x[2])
				max=a[i].x[2];
		}

		//prepare the graph
		mglGraph gr;
		gr.Rotate(50, 60);
		gr.SetRanges(min, max, min, max, min, max);
		//gr.SetOrigin(0, 0);

		//plot crystal BasisVectors
		int point1[12]={0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6};
		int point2[12]={1, 2, 3, 4, 5, 4, 6, 5, 6, 7, 7, 7};
		for(int i=0; i<12; i++)
		{
			mglData x(2), y(2), z(2);
			x.a[0]=a[point1[i]].x[0];
			y.a[0]=a[point1[i]].x[1];
			z.a[0]=a[point1[i]].x[2];
			x.a[1]=a[point2[i]].x[0];
			y.a[1]=a[point2[i]].x[1];
			z.a[1]=a[point2[i]].x[2];
			gr.Plot(x, y, z, "k-");
		}
		//plot particles
		mglData x(list->NumParticle()), y(list->NumParticle()), z(list->NumParticle());
		for(int i=0; i<list->NumParticle(); i++)
		{
			GeometryVector now=list->GetCartesianCoordinates(i);
			x.a[i]=now.x[0];
			y.a[i]=now.x[1];
			z.a[i]=now.x[2];
		}
		gr.Plot(x, y, z, "b o");

		prefix += std::string(FigureFormat);
		gr.WriteFrame(prefix.c_str());
		if(AlsoWriteEPS)
		{
			prefix.resize(prefix.size()-4);
			prefix+=std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else
	{
		std::cerr<<"error in Plots.cpp : Plot : unsupported dimension\n";
	}
#endif
}

const double NonLogScale_const = 20;
void PlotG1_2D(const std::vector<Configuration> & vc, long numBinsPerSide, std::string prefix, std::string remark)
{
#ifdef USE_MATHGL
	if (vc.size() == 0)
	{
		std::cerr << "Error in PlotG1_2D: no configuration!\n";
		return;
	}
	Configuration c = vc[0];
	GeometryVector a0(0, 0);
	GeometryVector a1 = a0 + c.GetBasisVector(0);
	GeometryVector a2 = a0 + c.GetBasisVector(1);
	GeometryVector a3 = a2 + c.GetBasisVector(0);

	//determine range of the plot
	double min = 0, max = 0;
	if (min > a0.x[0])
		min = a0.x[0];
	if (min > a0.x[1])
		min = a0.x[1];
	if (max < a0.x[0])
		max = a0.x[0];
	if (max < a0.x[1])
		max = a0.x[1];
	if (min > a1.x[0])
		min = a1.x[0];
	if (min > a1.x[1])
		min = a1.x[1];
	if (max < a1.x[0])
		max = a1.x[0];
	if (max < a1.x[1])
		max = a1.x[1];
	if (min > a2.x[0])
		min = a2.x[0];
	if (min > a2.x[1])
		min = a2.x[1];
	if (max < a2.x[0])
		max = a2.x[0];
	if (max < a2.x[1])
		max = a2.x[1];
	if (min > a3.x[0])
		min = a3.x[0];
	if (min > a3.x[1])
		min = a3.x[1];
	if (max < a3.x[0])
		max = a3.x[0];
	if (max < a3.x[1])
		max = a3.x[1];

	//prepare data
	mglData x(numBinsPerSide+1, numBinsPerSide+1);
	mglData y = x;
	mglData z = x;
	for(int i=0; i<numBinsPerSide+1; i++)
		for (int j = 0; j < numBinsPerSide + 1; j++)
		{
			GeometryVector temp = (double)(i) / numBinsPerSide*a1 + (double)(j) / numBinsPerSide*a2;
			x.SetVal(temp.x[0], i, j);
			y.SetVal(temp.x[1], i, j);
			z.SetVal(0.0, i, j);
		}

	size_t sumNumParticle = 0;
	for (auto iter = vc.begin(); iter != vc.end(); ++iter)
	{
		const Configuration & c = *iter;
		for (int i = 0; i < c.NumParticle(); i++)
		{
			GeometryVector t = c.GetRelativeCoordinates(i);
			int nx = std::floor(t.x[0] * numBinsPerSide);
			int ny = std::floor(t.x[1] * numBinsPerSide);
			z.SetVal(z.GetVal(nx, ny) + 1.0, nx, ny);
		}
		sumNumParticle += c.NumParticle();
	}
	for (int i = 0; i < numBinsPerSide; i++)
		z.SetVal(z.GetVal(i, 0), i, numBinsPerSide);
	for (int i = 0; i < numBinsPerSide + 1; i++)
		z.SetVal(z.GetVal(0, i), numBinsPerSide, i);
	for (int i = 0; i<numBinsPerSide + 1; i++)
		for (int j = 0; j < numBinsPerSide + 1; j++)
			z.SetVal(z.GetVal(i, j)/sumNumParticle*numBinsPerSide*numBinsPerSide, i, j);

	//generate ticks
	double maxS = z.Maximal();
	std::vector<double> ticks;
	std::stringstream ss;
	ticks.push_back(0);
	ticks.push_back(std::log10(0.1 * NonLogScale_const + 1));
	ticks.push_back(std::log10(0.2 * NonLogScale_const + 1));
	ticks.push_back(std::log10(0.5 * NonLogScale_const + 1));
	ticks.push_back(std::log10(1.0 * NonLogScale_const + 1));
	ss << "0\n0.1\n0.2\n0.5\n1";
	for (double i = 2; i<0.45*maxS; i *= 2)
	{
		ticks.push_back(std::log10(i * NonLogScale_const + 1));
		ss << '\n' << i;
	}
	ticks.push_back(std::log10(maxS * NonLogScale_const + 1));
	ss << '\n' << std::ceil(maxS);
	mglData Ticks(ticks.size());
	for (int i = 0; i<ticks.size(); i++)
		Ticks.SetVal(ticks[i], i);

	//non-linear transformation on z
	for (int i = 0; i<numBinsPerSide + 1; i++)
		for (int j = 0; j < numBinsPerSide + 1; j++)
			z.SetVal(std::log10(z.GetVal(i, j) * NonLogScale_const + 1), i, j);

	//prepare the graph
	mglGraph gr;

	gr.SetTickSkip(false);

	gr.SetSize(3200, 3200);
	//gr.SetFunc("","","","lg(c*10+1)");
	//gr.Title("S(k)");
	gr.Puts(mglPoint(0.0, 1.1), "g_1", ":C", 7);

	gr.SetRanges(min, max, min, max, 0, std::log10(std::ceil(maxS) * NonLogScale_const + 1));

	gr.Aspect(1.0, 1.0);
	gr.SetTicksVal('c', Ticks, ss.str().c_str());
	gr.Colorbar(">bcyr");

	gr.Puts(mglPoint(0.5*(max + min), min + 1.1*(max - min)), remark.c_str());

	//plot crystal BasisVectors
	{
		mglData x(2), y(2);
		x.a[0] = a0.x[0];
		y.a[0] = a0.x[1];
		x.a[1] = a1.x[0];
		y.a[1] = a1.x[1];
		gr.Plot(x, y, "k-");
		x.a[0] = a1.x[0];
		y.a[0] = a1.x[1];
		x.a[1] = a3.x[0];
		y.a[1] = a3.x[1];
		gr.Plot(x, y, "k-");
		x.a[0] = a2.x[0];
		y.a[0] = a2.x[1];
		x.a[1] = a3.x[0];
		y.a[1] = a3.x[1];
		gr.Plot(x, y, "k-");
		x.a[0] = a2.x[0];
		y.a[0] = a2.x[1];
		x.a[1] = a0.x[0];
		y.a[1] = a0.x[1];
		gr.Plot(x, y, "k-");
	}


	gr.Surf(x, y, z, "bcyr");

	//gr.Label('y', "y", 0);
	//gr.Label('x', "x", 0);

	prefix += std::string(FigureFormat);
	gr.WriteFrame(prefix.c_str());
	if (AlsoWriteEPS)
	{
		prefix.resize(prefix.size() - 4);
		prefix += std::string(".eps");
		gr.WriteEPS(prefix.c_str());
	}
#endif

}

#include "StructureFactor.h"
void PlotDirectionalStructureFactor2D(const std::vector<Configuration> & vc, long limit, std::string prefix, bool LogScale, double MinSOverride, double MaxSOverride, double TileSize)
{
	if (vc.empty())
	{
		std::cerr << "Error in PlotDirectionalStructureFactor2D: no configurations, exiting.\n";
		return;
	}

#ifdef USE_MATHGL
	GeometryVector v[2]={vc[0].GetReciprocalBasisVector(0), vc[0].GetReciprocalBasisVector(1)};
	double max=0;
	double maxS=0;
	double minS=0;

	//determine max and min of the graph
	for(long i=0; i<limit; i++)
	{
		for(long j=0; j<limit; j++)
		{
			GeometryVector k=static_cast<double>(i)*v[0]+static_cast<double>(j)*v[1];
			double tempMax = std::max(std::abs(k.x[0]), std::abs(k.x[1]));
			if(max<tempMax)
				max=tempMax;
		}
	}

	//determine the true limit to cover the canvas
	double truelimit = limit;
	PeriodicCellList<Empty> recip(2, v, max, false);
	recip.Insert(Empty(), GeometryVector(2));
	auto IterateFunction = [&truelimit, &max](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
	{
		if( std::abs(shift.x[0])<max && std::abs(shift.x[1])<max )
		{
			long l=std::abs(PeriodicShift[0]);
			if(truelimit<l)
				truelimit=l;
			l=std::abs(PeriodicShift[1]);
			if(truelimit<l)
				truelimit=l;
		}
	};
	recip.IterateThroughNeighbors(GeometryVector(2), std::sqrt(2.0)*max, IterateFunction);

	long siz=2*truelimit+1;
	mglData x(siz, siz);
	mglData y=x;
	mglData z=x;
	mglData r=x;
	GeometryVector k1=v[0]+v[1];
	for(long i=0; i<siz; i++)
	{
		for(long j=0; j<siz; j++)
		{
			GeometryVector k=static_cast<double>(i-truelimit)*v[0]+static_cast<double>(j-truelimit)*v[1];
			GeometryVector k2=k-0.5*k1;
			x.SetVal(k2.x[0], i, j);
			y.SetVal(k2.x[1], i, j);

			double s=0.0;
			for (auto iter = vc.begin(); iter != vc.end(); ++iter)
			{
				if (LogScale)
					s += std::log10(StructureFactor(*iter, k));
				else
					s += StructureFactor(*iter, k);
			}
			s /= vc.size();
			//if((!LogScale) && i-limit==0 && j-limit==0)
			//	s=std::sqrt(-1);//without log scale, plotting the foward scaterring will dramatically increase maxS
			if(i-truelimit==0 && j-truelimit==0)
				s=std::sqrt(-1);

			if (MinSOverride != 0 && s < MinSOverride)
				s = MinSOverride;
			if (MaxSOverride != 0 && s > MaxSOverride)
				s = MaxSOverride;

			//debug temp
			//std::cout<<s<<'\n';

			if (LogScale)
				z.SetVal(s, i, j);
			else
				z.SetVal(std::log10(NonLogScale_const * s + 1), i, j);

			if(maxS<s )
				maxS=s;
			if(minS>s)
				minS=s;
		}
	}
	if (MinSOverride != 0)
		minS = MinSOverride;
	if (MaxSOverride != 0)
		maxS = MaxSOverride;

	for (long i = 0; i<siz; i++)
		for(long j=0; j<siz; j++)
			r.SetVal(0.2*(maxS+minS)*TileSize, i, j);

	mglGraph gr;

	if(LogScale)
	{
		gr.SetSize(800, 600);
		gr.Title("log_{10}(S(k))");
		gr.SetTickRotate(false);
		gr.SetRanges( (-1.0)*max, max, (-1.0)*max, max, minS, std::ceil(maxS));
		gr.Aspect(0.75, 1.0);
		gr.Colorbar(">bcyr");

		if(limit<=50)
			gr.TileS(x, y, z, r, "bcyr");
		else
			gr.Surf(x, y, z, "bcyr");
		gr.Box();
		gr.Axis("U");

		gr.Label('y', "k_y", 0, "size 6");
		gr.Label('x', "k_x", 0, "value -0.1; size 6");
		//gr.Label('y', "log_{10}[S(k)]", 0, "value -3.3; size 5");
	}
	else
	{
		//generate ticks
		std::vector<double> ticks;
		std::stringstream ss;
		ticks.push_back(0);
		ticks.push_back(std::log10(0.1 * NonLogScale_const + 1));
		ticks.push_back(std::log10(0.2 * NonLogScale_const + 1));
		ticks.push_back(std::log10(0.5 * NonLogScale_const + 1));
		ticks.push_back(std::log10(1.0 * NonLogScale_const + 1));
		ss<<"0\n0.1\n0.2\n0.5\n1";
		for(double i=2; i<0.45*maxS; i*=2)
		{
			ticks.push_back(std::log10(i * NonLogScale_const + 1));
			ss<<'\n'<<i;
		}
		ticks.push_back(std::log10(maxS * NonLogScale_const + 1));
		ss<<'\n'<<std::ceil(maxS);
		//ticks.push_back(0.0);
		//ss << "0\n";
		//ticks.push_back(maxS);
		//ss << "\\infty\n";

		mglData Ticks(ticks.size());
		for(int i=0; i<ticks.size(); i++)
			Ticks.SetVal(ticks[i], i);

		gr.SetTickSkip(false);

		gr.SetSize(800, 800);
		//gr.SetFunc("","","","lg(c*10+1)");
		//gr.Title("S(k)");
		gr.Puts(mglPoint(0.0, 1.1), "S(k)", ":C", 7);

		gr.SetRanges((-1.0)*max, max, (-1.0)*max, max, 0, std::log10(std::ceil(maxS) * NonLogScale_const + 1));

		gr.Aspect(1.0, 1.0);
		gr.Box();
		gr.SetTicksVal('c', Ticks, ss.str().c_str());
		gr.Colorbar(">bcyr");
		gr.Axis();

		if (limit <= 50)
			gr.TileS(x, y, z, r, "bcyr");
		else
			gr.Surf(x, y, z, "bcyr");

		gr.Label('y', "k_2", 0);
		gr.Label('x', "k_1", 0);
	}

	prefix += std::string(FigureFormat);
	gr.WriteFrame(prefix.c_str());
	if(AlsoWriteEPS)
	{
		prefix.resize(prefix.size()-4);
		prefix+=std::string(".eps");
		gr.WriteEPS(prefix.c_str());
	}
#endif
}
void PlotDirectionalStructureFactor2D(const Configuration & c, long limit, std::string prefix, bool LogScale, double MinSOverride, double MaxSOverride, double TileSize)
{
	PlotDirectionalStructureFactor2D(std::vector<Configuration>(1, c), limit, prefix, LogScale, MinSOverride, MaxSOverride, TileSize);
}

void Plots(const ParameterSet & Param, const std::vector<LatticeSumSolver *> * Solvers, const char * Prefix, size_t Enphasis)
{
	::IsotropicPairPotential * pPot = Param.GetPotential();

	std::string filename(Prefix);
	filename+="_Parameters.txt";
	std::fstream ofile;
	ofile.open(filename.c_str(), std::fstream::out);
	ofile.precision(16);
	Param.Write(ofile);
	ofile.close();

	::Plots(pPot, Solvers, Prefix, Enphasis);

	delete pPot;
}

void Plots(IsotropicPairPotential * pPot, const std::vector<LatticeSumSolver *> * Solvers, const char * Prefix, size_t Enphasis)
{
	std::fstream ofile;

	std::string filename(Prefix);

	filename.assign(Prefix);
	filename+="_Potential.txt";
	ofile.open(filename.c_str(), std::fstream::out);
#ifdef USE_MATHGL
	std::vector<double> potentialx, potentialy;
#endif
	double LowestPotential=::MaxEnergy;
	double HighestPotential = 0;
	bool startTrackHighest = false;
	double LowestNumber=0;
	double ShowMinDist=0.5;
	if(Solvers!=nullptr)
		if(Solvers->size()!=0)
		{
			double tempDist=1;
			while(Solvers->at(0)->Terms.size()<2)
			{
				Solvers->at(0)->UpdateTerms(tempDist);
				tempDist*=2;
			}
			ShowMinDist=Solvers->at(0)->Terms[1].distance;
		}
		for(double i=0.001; i<pPot->Rcut; i+=0.001)
		{
			double result=pPot->IsotropicPairEnergy(i, "", "");
			ofile<<i<<" \t"<<result<<'\n';
			if(result<LowestPotential)
			{
				LowestPotential=result;
#ifdef USE_MATHGL
				LowestNumber=potentialx.size();
#endif
			}
			if(result>LowestPotential && LowestPotential<MaxEnergy)
				startTrackHighest=true;
			if(i>ShowMinDist)
				startTrackHighest=true;
			if(startTrackHighest)
			{
				if(result>HighestPotential)
					HighestPotential=result;
			}
#ifdef USE_MATHGL
			potentialx.push_back(i);
			potentialy.push_back(result);
#endif
		}
		HighestPotential*=2.001;

		ofile.close();
#ifdef USE_MATHGL
		{
			mglData x, y;
			x.Set(&potentialx[0], potentialx.size());
			y.Set(&potentialy[0], potentialy.size());
			mglGraph gr;
			gr.SetRanges(0, pPot->Rcut, LowestPotential, HighestPotential);
			gr.SetOrigin(0, LowestPotential);
			gr.SetFontSize(6);
			gr.SetTickTempl('y', " ");
			gr.SetTickShift(mglPoint(0.1, 0, 0));
			gr.Axis();
			gr.SetTickShift(mglPoint(0.0, 0.0));
			gr.SetTickTempl('y', "");
			gr.SetTickTempl('x', " ");
			gr.Axis();
			gr.Label('y', "u_2(r)", 0);
			gr.Label('x', "r", 0);
			gr.SetRanges(0, pPot->Rcut, LowestPotential, HighestPotential);
			gr.SetOrigin(pPot->Rcut, HighestPotential);
			gr.SetTickTempl('x', " ");
			gr.SetTickTempl('y', " ");
			gr.Axis();
			gr.Plot(x, y);
			filename.assign(Prefix);
			filename+="_Potential";
			filename += FigureFormat;
			gr.WriteFrame(filename.c_str());
			if(AlsoWriteEPS)
			{
				filename.assign(Prefix);
				filename+="_Potential.eps";
				gr.WriteEPS(filename.c_str());
			}
		}
#endif

		if(Solvers!=nullptr)
		{
			filename.assign(Prefix);
			filename+="_LatticeSum.txt";
			std::vector<const char *> tags;
			std::vector<double> volumes;
			std::vector<std::vector<double>> results;
			double LowestSum=::MaxEnergy;
			auto temp=(*Solvers)[0]->GetStructure();
			double HighestSum=3*std::abs(pPot->Energy(temp));
			//calculation
			for(auto iter=Solvers->begin(); iter!=Solvers->end(); iter++)
			{
				tags.push_back((*iter)->Tag());
				results.push_back(std::vector<double>());

				bool startTrackHighest=false;

				for(double Volume=0.3; Volume<10; Volume /=0.998)
				{
					if(iter==Solvers->begin())
						volumes.push_back(Volume);
					double sum=(*iter)->LatticeSum(Volume, *pPot);
					if(sum<LowestSum)
					{
						LowestSum=sum;
					}
					results.back().push_back(sum);
				}

			}


			//output
			ofile.open(filename, std::fstream::out);
			ofile.precision(10);
			ofile<<std::scientific;
			ofile<<"1/Density\t";
			for(size_t i=0; i<tags.size(); i++)
			{
				ofile<<tags[i]<<" \t";
			}
			ofile<<'\n';
			for(size_t i=0; i<volumes.size(); i++)
			{
				ofile<<volumes[i]<<" \t";
				for(size_t j=0; j<tags.size(); j++)
					ofile<<results[j][i]<<" \t";
				ofile<<'\n';
			}

			ofile.close();
#ifdef USE_MATHGL
			{
				mglGraph gr;
				gr.SetRanges(0, pPot->Rcut, LowestSum, HighestSum);
				gr.Axis();
				gr.Label('y', "E/N", 0);
				gr.Label('x', "1/Density", 0);

				mglData x, y;
				x.Set(&volumes[0], volumes.size());

				y.Set(&results[0][0], results[0].size());
				gr.Plot(x, y, "k1-3");

				for(size_t j=1; j<tags.size(); j++)
				{
					y.Set(&results[j][0], results[j].size());
					if(j!=Enphasis)
						gr.Plot(x, y, "C8-1");
					else
						gr.Plot(x, y, "r8-2");
				}
				filename.assign(Prefix);
				filename+="_LatticeSum";
				filename += FigureFormat;
				gr.WriteFrame(filename.c_str());
				if(AlsoWriteEPS)
				{
					filename.assign(Prefix);
					filename+="_LatticeSum.eps";
					gr.WriteEPS(filename.c_str());
				}
			}
#endif

		}

}


void PlotPath(std::string prefix, const Configuration & List, const std::string & remark, const std::vector<size_t> & path)
{
	const Configuration * list= & List;
#ifdef USE_MATHGL
	GeometryVector basis0=list->GetBasisVector(0);
	GeometryVector basis1=list->GetBasisVector(1);

	if(list->GetDimension()==2)
	{
		//calculate the range of the plot
		GeometryVector a0(0, 0);
		GeometryVector a1=a0+basis0;
		GeometryVector a2=a0+basis1;
		GeometryVector a3=a2+basis0;
		double min=0, max=0;
		if(min>a0.x[0])
			min=a0.x[0];
		if(min>a0.x[1])
			min=a0.x[1];
		if(max<a0.x[0])
			max=a0.x[0];
		if(max<a0.x[1])
			max=a0.x[1];
		if(min>a1.x[0])
			min=a1.x[0];
		if(min>a1.x[1])
			min=a1.x[1];
		if(max<a1.x[0])
			max=a1.x[0];
		if(max<a1.x[1])
			max=a1.x[1];
		if(min>a2.x[0])
			min=a2.x[0];
		if(min>a2.x[1])
			min=a2.x[1];
		if(max<a2.x[0])
			max=a2.x[0];
		if(max<a2.x[1])
			max=a2.x[1];
		if(min>a3.x[0])
			min=a3.x[0];
		if(min>a3.x[1])
			min=a3.x[1];
		if(max<a3.x[0])
			max=a3.x[0];
		if(max<a3.x[1])
			max=a3.x[1];

		//prepare the graph
		mglGraph gr;
		gr.SetRanges(min, max, min, max);
		gr.Aspect(1.0, 1.0);
		gr.SetSize(800, 800);
		//gr.SetOrigin(0, 0);
		gr.Puts(mglPoint(0.5*(max+min), min+1.1*(max-min)), remark.c_str());


		//plot crystal BasisVectors
		{
			mglData x(2), y(2);
			x.a[0]=a0.x[0];
			y.a[0]=a0.x[1];
			x.a[1]=a1.x[0];
			y.a[1]=a1.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a1.x[0];
			y.a[0]=a1.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a0.x[0];
			y.a[1]=a0.x[1];
			gr.Plot(x, y, "k-");
		}
		//plot particles
		double radius=0.1*(max-min)/std::sqrt(list->NumParticle());
		for(int i=0; i<list->NumParticle(); i++)
		{
			GeometryVector now=list->GetCartesianCoordinates(i);
			PlotDisk2D(gr, now.x[0], now.x[1], radius, "b", 10);
		}

		//plot path
		for(auto iter=path.begin()+1; iter<path.end(); iter++)
		{
			mglData x(2), y(2);
			GeometryVector v1=list->GetCartesianCoordinates(*(iter-1));
			x[0]=v1.x[0];
			y[0]=v1.x[1];
			GeometryVector v2=list->GetCartesianCoordinates(*iter);
			GeometryVector delta=list->GetRelativeCoordinates(*iter)-list->GetRelativeCoordinates(*(iter-1));
			double f[2] = { std::floor(delta.x[0] + 0.5), std::floor(delta.x[1] + 0.5) };
			delta.x[0]-=f[0];
			delta.x[1]-=f[1];
			
			delta = list->RelativeCoord2CartesianCoord(delta);

			v2=delta+v1;
			x[1]=v2.x[0];
			y[1]=v2.x[1];
			gr.Plot(x, y, "k");

			if (f[0] != 0.0 || f[1] != 0.0)
			{
				v2 = list->GetCartesianCoordinates(*iter);
				v1 = v2 - delta;
				x[0] = v1.x[0];
				y[0] = v1.x[1];
				x[1] = v2.x[0];
				y[1] = v2.x[1];
				gr.Plot(x, y, "k");
			}
		}

		prefix += std::string(FigureFormat);
		gr.WriteFrame(prefix.c_str());
		if(AlsoWriteEPS)
		{
			prefix.resize(prefix.size()-4);
			prefix+=std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else
	{
		std::cerr<<"error in Plots.cpp : Plot : unsupported dimension\n";
	}
#endif
}


//given a configuration and a "deformation" (i.e. a displacement on each particle), plot it
void PlotDisplacement(std::string prefix, const Configuration & List, const std::string & remark, const std::vector<double> & displacement)
{
	const Configuration * list= & List;
#ifdef USE_MATHGL
	GeometryVector basis0=list->GetBasisVector(0);
	GeometryVector basis1=list->GetBasisVector(1);

	if(list->GetDimension()==2)
	{
		//calculate the range of the plot
		GeometryVector a0(0, 0);
		GeometryVector a1=a0+basis0;
		GeometryVector a2=a0+basis1;
		GeometryVector a3=a2+basis0;
		double min=0, max=0;
		if(min>a0.x[0])
			min=a0.x[0];
		if(min>a0.x[1])
			min=a0.x[1];
		if(max<a0.x[0])
			max=a0.x[0];
		if(max<a0.x[1])
			max=a0.x[1];
		if(min>a1.x[0])
			min=a1.x[0];
		if(min>a1.x[1])
			min=a1.x[1];
		if(max<a1.x[0])
			max=a1.x[0];
		if(max<a1.x[1])
			max=a1.x[1];
		if(min>a2.x[0])
			min=a2.x[0];
		if(min>a2.x[1])
			min=a2.x[1];
		if(max<a2.x[0])
			max=a2.x[0];
		if(max<a2.x[1])
			max=a2.x[1];
		if(min>a3.x[0])
			min=a3.x[0];
		if(min>a3.x[1])
			min=a3.x[1];
		if(max<a3.x[0])
			max=a3.x[0];
		if(max<a3.x[1])
			max=a3.x[1];

		//prepare the graph
		mglGraph gr;
		gr.SetRanges(min, max, min, max);
		gr.Aspect(1.0, 1.0);
		gr.SetSize(3200, 3200);
		//gr.SetOrigin(0, 0);
		gr.Puts(mglPoint(0.5*(max+min), min+1.1*(max-min)), remark.c_str());


		//plot crystal BasisVectors
		{
			mglData x(2), y(2);
			x.a[0]=a0.x[0];
			y.a[0]=a0.x[1];
			x.a[1]=a1.x[0];
			y.a[1]=a1.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a1.x[0];
			y.a[0]=a1.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a0.x[0];
			y.a[1]=a0.x[1];
			gr.Plot(x, y, "k-");
		}
		//plot particles
		//mglData x(list->NumParticle()), y(list->NumParticle());
		double radius=0.1*(max-min)/std::sqrt(list->NumParticle());
		for(int i=0; i<list->NumParticle(); i++)
		{
			GeometryVector now=list->GetCartesianCoordinates(i);
			//x.a[i]=now.x[0];
			//y.a[i]=now.x[1];
			PlotDisk2D(gr, now.x[0], now.x[1], radius, "b", 10);

			//if(list->GetCharacteristics(i)==AtomInfo("A"))
			//	PlotDisk2D(gr, now.x[0], now.x[1], 0.5, "b", 10);
			//else
			//	PlotDisk2D(gr, now.x[0], now.x[1], 0.7, "r", 10);
		}
		//gr.Plot(x, y, "b 8.");

		//determine a rescaling of lines
		double TypicalLength = std::pow( list->PeriodicVolume()/list->NumParticle(), (1.0)/list->GetDimension() );
		double LargestLength2 = 0.0;
		for(int i=0; i<list->NumParticle(); i++)
		{
			GeometryVector v(displacement[2*i], displacement[2*i+1]);
			if(v.Modulus2()>LargestLength2)
				LargestLength2=v.Modulus2();
		}
		//double coeff = TypicalLength/std::sqrt(LargestLength2)*0.8;
		double coeff = TypicalLength/std::sqrt(LargestLength2)*3.0;
		//plot path
		gr.SetArrowSize(0.2);
		for(int i=0; i<list->NumParticle(); i++)
		{
			GeometryVector v1=list->GetCartesianCoordinates(i);
			//mglData x(2), y(2);
			//x[0]=v1.x[0];
			//y[0]=v1.x[1];
			//x[1]=x[0]+coeff*displacement[2*i];
			//y[1]=y[0]+coeff*displacement[2*i+1];
			//gr.Plot(x, y, "k");
			mglPoint p1, p2;
			p1.x=v1.x[0];
			p1.y=v1.x[1];
			p2.x=v1.x[0]+coeff*displacement[2*i];
			p2.y=v1.x[1]+coeff*displacement[2*i+1];
			gr.Line(p1, p2, "kA");
		}
		//temp
		//mglPoint p1, p2;
		//p1.x=0.2;
		//p1.y=0.2;
		//p2.x=0.2;
		//p2.y=0.8;
		//gr.Line(p1, p2, "B");
		//p1.x=0.4;
		//p1.y=0.2;
		//p2.x=0.4;
		//p2.y=0.8;
		//gr.Line(p1, p2, "B");
		//p1.x=0.6;
		//p1.y=0.2;
		//p2.x=0.6;
		//p2.y=0.8;
		//gr.Line(p1, p2, "B");
		//p1.x=0.8;
		//p1.y=0.2;
		//p2.x=0.8;
		//p2.y=0.8;
		//gr.Line(p1, p2, "B");
		//p1.x=0.2;
		//p1.y=0.8;
		//p2.x=0.8;
		//p2.y=0.8;
		//gr.Line(p1, p2, "B");
		//p1.x=0.2;
		//p1.y=0.2;
		//p2.x=0.8;
		//p2.y=0.2;
		//gr.Line(p1, p2, "B");
		//p1.x=0.2;
		//p1.y=0.2;
		//p2.x=0.8;
		//p2.y=0.2;
		//gr.Line(p1, p2, "B");
		//p1.x=0.2;
		//p1.y=0.4;
		//p2.x=0.8;
		//p2.y=0.4;
		//gr.Line(p1, p2, "B");
		//p1.x=0.2;
		//p1.y=0.6;
		//p2.x=0.8;
		//p2.y=0.6;
		//gr.Line(p1, p2, "B");

		prefix += std::string(FigureFormat);
		gr.WriteFrame(prefix.c_str());
		if(AlsoWriteEPS)
		{
			prefix.resize(prefix.size()-4);
			prefix+=std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else
	{
		std::cerr<<"error in Plots.cpp : Plot : unsupported dimension\n";
	}
#endif
}


void PlotPacking2D(std::string prefix, const SpherePacking & List, const std::string & remark, std::function<std::string(size_t Num)> GetMGLOptionsFunc)
{
	const SpherePacking * list= & List;
#ifdef USE_MATHGL

	if(list->GetDimension()==2)
	{
		//calculate the range of the plot
		GeometryVector a0(0, 0);
		GeometryVector a1=a0+list->GetBasisVector(0);
		GeometryVector a2=a0+list->GetBasisVector(1);
		GeometryVector a3=a2+list->GetBasisVector(0);
		double min=0, max=0;
		if(min>a0.x[0])
			min=a0.x[0];
		if(min>a0.x[1])
			min=a0.x[1];
		if(max<a0.x[0])
			max=a0.x[0];
		if(max<a0.x[1])
			max=a0.x[1];
		if(min>a1.x[0])
			min=a1.x[0];
		if(min>a1.x[1])
			min=a1.x[1];
		if(max<a1.x[0])
			max=a1.x[0];
		if(max<a1.x[1])
			max=a1.x[1];
		if(min>a2.x[0])
			min=a2.x[0];
		if(min>a2.x[1])
			min=a2.x[1];
		if(max<a2.x[0])
			max=a2.x[0];
		if(max<a2.x[1])
			max=a2.x[1];
		if(min>a3.x[0])
			min=a3.x[0];
		if(min>a3.x[1])
			min=a3.x[1];
		if(max<a3.x[0])
			max=a3.x[0];
		if(max<a3.x[1])
			max=a3.x[1];

		//prepare the graph
		mglGraph gr;
		gr.SetRanges(min, max, min, max);
		gr.Aspect(1.0, 1.0);
		gr.SetSize(800, 800);
		//gr.SetOrigin(0, 0);
		gr.Puts(mglPoint(0.5*(max+min), min+1.1*(max-min)), remark.c_str());


		//plot crystal BasisVectors
		{
			mglData x(2), y(2);
			x.a[0]=a0.x[0];
			y.a[0]=a0.x[1];
			x.a[1]=a1.x[0];
			y.a[1]=a1.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a1.x[0];
			y.a[0]=a1.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a0.x[0];
			y.a[1]=a0.x[1];
			gr.Plot(x, y, "k-");
		}
		//plot particles
		for(int i=0; i<list->NumParticle(); i++)
		{
			double r=list->GetCharacteristics(i);
			GeometryVector now=list->GetCartesianCoordinates(i);
			::PlotDisk2D(gr, now.x[0], now.x[1], r, GetMGLOptionsFunc(i).c_str(), 50);

			//images
			//{
			//	GeometryVector n2=now+list->GetBasisVector(0);
			//	::PlotDisk2D(gr, n2.x[0], n2.x[1], r, "b");
			//}
			//{
			//	GeometryVector n2=now+list->GetBasisVector(1);
			//	::PlotDisk2D(gr, n2.x[0], n2.x[1], r, "b");
			//}
			//{
			//	GeometryVector n2=now-list->GetBasisVector(0);
			//	::PlotDisk2D(gr, n2.x[0], n2.x[1], r, "b");
			//}
			//{
			//	GeometryVector n2=now-list->GetBasisVector(1);
			//	::PlotDisk2D(gr, n2.x[0], n2.x[1], r, "b");
			//}
			//{
			//	GeometryVector n2=now+list->GetBasisVector(0)+list->GetBasisVector(1);
			//	::PlotDisk2D(gr, n2.x[0], n2.x[1], r, "b");
			//}
			//{
			//	GeometryVector n2=now-list->GetBasisVector(0)+list->GetBasisVector(1);
			//	::PlotDisk2D(gr, n2.x[0], n2.x[1], r, "b");
			//}
			//{
			//	GeometryVector n2=now-list->GetBasisVector(0)-list->GetBasisVector(1);
			//	::PlotDisk2D(gr, n2.x[0], n2.x[1], r, "b");
			//}
			//{
			//	GeometryVector n2=now+list->GetBasisVector(0)-list->GetBasisVector(1);
			//	::PlotDisk2D(gr, n2.x[0], n2.x[1], r, "b");
			//}
		}

		prefix += std::string(FigureFormat);
		gr.WriteFrame(prefix.c_str());
		if(AlsoWriteEPS)
		{
			prefix.resize(prefix.size()-4);
			prefix+=std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else
	{
		std::cerr<<"error in Plots.cpp : PlotPlotPacking2D : unsupported dimension\n";
	}
#endif
}

template<typename add>
void PlotPacking2D(std::string prefix, const PeriodicCellList<add> & List, const std::string & remark, double r)
{
	PlotPacking2D(prefix, SpherePacking(List, r), remark);
}


void gracePlot::autoScale()
{
	if (dataSets.size() == 0)
	{
		std::cerr << "Error in gracePlot : no data set, cannot auto scale!\n";
		return;
	}
	MinX = dataSets[0].result[0].x[0], MaxX = dataSets[0].result[0].x[0];
	MinY = dataSets[0].result[0].x[1], MaxY = dataSets[0].result[0].x[1];
	for (dataSet & set : dataSets)
	{
		const std::vector<GeometryVector> & result = set.result;
		const size_t NumData = result.size();
		for (size_t i = 0; i<NumData; i++)
		{
			if (MinX>result[i].x[0])
				MinX = result[i].x[0];
			if (MaxX<result[i].x[0])
				MaxX = result[i].x[0];
			if (MinY>result[i].x[1])
				MinY = result[i].x[1];
			if (MaxY<result[i].x[1])
				MaxY = result[i].x[1];
		}
	}
	if (xScale == "Logarithmic" && MinX <= 0.0)
		MinX = 1e-6*MaxX;
	if (yScale == "Logarithmic" && MinY <= 0.0)
		MinY = 1e-6*MaxY;
}
void gracePlot::autoTick()
{
	if (xScale == "Logarithmic")
		TickX = 10;
	else
	{
		double TickX = std::pow(10.0, std::floor(std::log10(MaxX - MinX) - 0.1761));
		if ((MaxX - MinX) / TickX < 3)
			TickX *= 0.5;
		else if ((MaxX - MinX) / TickX>6)
			TickX *= 2;
	}

	if (yScale == "Logarithmic")
		TickY = 10;
	else
	{
		double TickY = std::pow(10.0, std::floor(std::log10(MaxY - MinY) - 0.1761));
		if ((MaxY - MinY) / TickY < 3)
			TickY *= 0.5;
		else if ((MaxY - MinY) / TickY>6)
			TickY *= 2;
	}
}
void gracePlot::autoScaleAndTick()
{
	autoScale();
	autoTick();
}
void gracePlot::addDataSet(std::vector<GeometryVector> data, std::string legend, int lineColor, int lineStyle, int symbolColor, int symbolStyle, std::string type, double symbolSize) // for colors and styles, -1 means auto, and 0 means none
{
	if (data.size() == 0)
	{
		std::cerr << "Error in gracePlot : no data! Exiting\n";
		return;
	}
	dataSet temp;
	std::swap(data, temp.result);
	temp.type = type;
	temp.legend = legend;

	if (lineColor == -1)
		temp.lineColor = (dataSets.size() % NColor) + 1;
	else if (lineColor == 0)
		temp.lineColor = 0;
	else
		temp.lineColor = ((lineColor-1) % NColor) + 1;

	if (lineStyle == -1)
		temp.lineStyle = (dataSets.size() % NLineStyle) + 1;
	else if (NLineStyle == 0)
		temp.lineStyle = 0;
	else
		temp.lineStyle = ((lineStyle-1) % NLineStyle) + 1;

	if (symbolColor == -1)
		temp.symbolColor = (dataSets.size() % NColor) + 1;
	else if (symbolColor == 0)
		temp.symbolColor = 0;
	else
		temp.symbolColor = ((symbolColor-1) % NColor) + 1;

	if (symbolStyle == -1 || symbolStyle == 0)
		temp.symbolStyle = 0;
	else
		temp.symbolStyle = ((symbolStyle-1) % NSymbol) + 1;
	temp.symbolSize = symbolSize;

	this->dataSets.push_back(temp);
}

void gracePlot::outputFigure(const std::string & OutFilePrefix)
{
	if (dataSets.size() == 0)
	{
		std::cerr << "Error in gracePlot : no data set! Exiting\n";
		return;
	}
	std::string name = OutFilePrefix + std::string(".agr");
	std::fstream ofile(name.c_str(), std::fstream::out);
	ofile.precision(15);
	ofile << "# Grace project file\n";
	ofile << "# Generated by gracePlot::outputFigure. For more information, contact Ge Zhang (gezhang@alumni.princeton.edu).\n";
	ofile << "@version 50122\n";
	ofile << "@page size 792, 612\n";
	ofile << "@page scroll 5%\n";
	ofile << "@page inout 5%\n";
	ofile << "@link page off\n";
	ofile << "@map font 0 to \"Times-Roman\", \"Times-Roman\"\n";
	ofile << "@map font 1 to \"Times-Italic\", \"Times-Italic\"\n";
	ofile << "@map font 2 to \"Times-Bold\", \"Times-Bold\"\n";
	ofile << "@map font 3 to \"Times-BoldItalic\", \"Times-BoldItalic\"\n";
	ofile << "@map font 4 to \"Helvetica\", \"Helvetica\"\n";
	ofile << "@map font 5 to \"Helvetica-Oblique\", \"Helvetica-Oblique\"\n";
	ofile << "@map font 6 to \"Helvetica-Bold\", \"Helvetica-Bold\"\n";
	ofile << "@map font 7 to \"Helvetica-BoldOblique\", \"Helvetica-BoldOblique\"\n";
	ofile << "@map font 8 to \"Courier\", \"Courier\"\n";
	ofile << "@map font 9 to \"Courier-Oblique\", \"Courier-Oblique\"\n";
	ofile << "@map font 10 to \"Courier-Bold\", \"Courier-Bold\"\n";
	ofile << "@map font 11 to \"Courier-BoldOblique\", \"Courier-BoldOblique\"\n";
	ofile << "@map font 12 to \"Symbol\", \"Symbol\"\n";
	ofile << "@map font 13 to \"ZapfDingbats\", \"ZapfDingbats\"\n";
	ofile << "@map color 0 to (255, 255, 255), \"white\"\n";
	ofile << "@map color 1 to (0, 0, 0), \"black\"\n";
	ofile << "@map color 2 to (255, 0, 0), \"red\"\n";
	ofile << "@map color 3 to (0, 255, 0), \"green\"\n";
	ofile << "@map color 4 to (0, 0, 255), \"blue\"\n";
	ofile << "@map color 5 to (188, 143, 143), \"brown\"\n";
	ofile << "@map color 6 to (220, 220, 220), \"grey\"\n";
	ofile << "@map color 7 to (148, 0, 211), \"violet\"\n";
	ofile << "@map color 8 to (0, 255, 255), \"cyan\"\n";
	ofile << "@map color 9 to (255, 0, 255), \"magenta\"\n";
	ofile << "@map color 10 to (255, 165, 0), \"orange\"\n";
	ofile << "@map color 11 to (114, 33, 188), \"indigo\"\n";
	ofile << "@map color 12 to (103, 7, 72), \"maroon\"\n";
	ofile << "@map color 13 to (64, 224, 208), \"turquoise\"\n";
	ofile << "@map color 14 to (0, 139, 0), \"green4\"\n";
	ofile << "@map color 15 to (255, 255, 0), \"yellow\"\n";
	ofile << "@reference date 0\n";
	ofile << "@date wrap off\n";
	ofile << "@date wrap year 1950\n";
	ofile << "@default linewidth 1.0\n";
	ofile << "@default linestyle 1\n";
	ofile << "@default color 1\n";
	ofile << "@default pattern 1\n";
	ofile << "@default font 0\n";
	ofile << "@default char size 1.000000\n";
	ofile << "@default symbol size 1.000000\n";
	ofile << "@default sformat \"%.8g\"\n";
	ofile << "@background color 0\n";
	ofile << "@page background fill off\n";
	ofile << "@timestamp off\n";
	ofile << "@timestamp 0.03, 0.03\n";
	ofile << "@timestamp color 1\n";
	ofile << "@timestamp rot 0\n";
	ofile << "@timestamp font 0\n";
	ofile << "@timestamp char size 1.000000\n";
	ofile << "@timestamp def \"";
	struct tm * now = localtime(&::ProgramStart);
	ofile << (now->tm_year + 1900) << '-'
		<< (now->tm_mon + 1) << '-'
		<< now->tm_mday << ' '
		<< now->tm_hour << ':'
		<< now->tm_min << ':'
		<< now->tm_sec;

	ofile << "\"\n";
	ofile << "@r0 off\n";
	ofile << "@link r0 to g0\n";
	ofile << "@r0 type above\n";
	ofile << "@r0 linestyle 1\n";
	ofile << "@r0 linewidth 1.0\n";
	ofile << "@r0 color 1\n";
	ofile << "@r0 line 0, 0, 0, 0\n";
	ofile << "@r1 off\n";
	ofile << "@link r1 to g0\n";
	ofile << "@r1 type above\n";
	ofile << "@r1 linestyle 1\n";
	ofile << "@r1 linewidth 1.0\n";
	ofile << "@r1 color 1\n";
	ofile << "@r1 line 0, 0, 0, 0\n";
	ofile << "@r2 off\n";
	ofile << "@link r2 to g0\n";
	ofile << "@r2 type above\n";
	ofile << "@r2 linestyle 1\n";
	ofile << "@r2 linewidth 1.0\n";
	ofile << "@r2 color 1\n";
	ofile << "@r2 line 0, 0, 0, 0\n";
	ofile << "@r3 off\n";
	ofile << "@link r3 to g0\n";
	ofile << "@r3 type above\n";
	ofile << "@r3 linestyle 1\n";
	ofile << "@r3 linewidth 1.0\n";
	ofile << "@r3 color 1\n";
	ofile << "@r3 line 0, 0, 0, 0\n";
	ofile << "@r4 off\n";
	ofile << "@link r4 to g0\n";
	ofile << "@r4 type above\n";
	ofile << "@r4 linestyle 1\n";
	ofile << "@r4 linewidth 1.0\n";
	ofile << "@r4 color 1\n";
	ofile << "@r4 line 0, 0, 0, 0\n";
	ofile << "@g0 on\n";
	ofile << "@g0 hidden false\n";
	ofile << "@g0 type XY\n";
	ofile << "@g0 stacked false\n";
	ofile << "@g0 bar hgap 0.000000\n";
	ofile << "@g0 fixedpoint off\n";
	ofile << "@g0 fixedpoint type 0\n";
	ofile << "@g0 fixedpoint xy 0.000000, 0.000000\n";
	ofile << "@g0 fixedpoint format general general\n";
	ofile << "@g0 fixedpoint prec 6, 6\n";
	ofile << "@with g0\n";
	ofile << "@    world " << MinX << ", " << MinY << ", " << MaxX << ", " << MaxY << '\n';
	ofile << "@    stack world 0, 0, 0, 0\n";
	ofile << "@    znorm 1\n";
	ofile << "@    view 0.250000, 0.200000, 1.150000, 0.850000\n";
	ofile << "@    title \"" << Title << "\"\n";
	ofile << "@    title font 0\n";
	ofile << "@    title size 1.500000\n";
	ofile << "@    title color 1\n";
	ofile << "@    subtitle \"\"\n";
	ofile << "@    subtitle font 0\n";
	ofile << "@    subtitle size 1.000000\n";
	ofile << "@    subtitle color 1\n";
	ofile << "@    xaxes scale " << xScale << "\n";
	ofile << "@    yaxes scale " << yScale << "\n";
	ofile << "@    xaxes invert off\n";
	ofile << "@    yaxes invert off\n";
	ofile << "@    xaxis  on\n";
	ofile << "@    xaxis  type zero false\n";
	ofile << "@    xaxis  offset 0.000000 , 0.000000\n";
	ofile << "@    xaxis  bar on\n";
	ofile << "@    xaxis  bar color 1\n";
	ofile << "@    xaxis  bar linestyle 1\n";
	ofile << "@    xaxis  bar linewidth 1.0\n";
	ofile << "@    xaxis  label \"" << xLabel << "\"\n";
	ofile << "@    xaxis  label layout para\n";
	ofile << "@    xaxis  label place auto\n";
	ofile << "@    xaxis  label char size 3.000000\n";
	ofile << "@    xaxis  label font 0\n";
	ofile << "@    xaxis  label color 1\n";
	ofile << "@    xaxis  label place normal\n";
	ofile << "@    xaxis  tick on\n";
	ofile << "@    xaxis  tick major " << TickX << "\n";
	ofile << "@    xaxis  tick minor ticks 0\n";
	ofile << "@    xaxis  tick default 6\n";
	ofile << "@    xaxis  tick place rounded true\n";
	ofile << "@    xaxis  tick in\n";
	ofile << "@    xaxis  tick major size 1.000000\n";
	ofile << "@    xaxis  tick major color 1\n";
	ofile << "@    xaxis  tick major linewidth 1.0\n";
	ofile << "@    xaxis  tick major linestyle 1\n";
	ofile << "@    xaxis  tick major grid off\n";
	ofile << "@    xaxis  tick minor color 1\n";
	ofile << "@    xaxis  tick minor linewidth 1.0\n";
	ofile << "@    xaxis  tick minor linestyle 1\n";
	ofile << "@    xaxis  tick minor grid off\n";
	ofile << "@    xaxis  tick minor size 0.500000\n";
	ofile << "@    xaxis  ticklabel on\n";
	ofile << "@    xaxis  ticklabel format general\n";
	ofile << "@    xaxis  ticklabel prec 5\n";
	ofile << "@    xaxis  ticklabel formula \"\"\n";
	ofile << "@    xaxis  ticklabel append \"\"\n";
	ofile << "@    xaxis  ticklabel prepend \"\"\n";
	ofile << "@    xaxis  ticklabel angle 0\n";
	ofile << "@    xaxis  ticklabel skip 0\n";
	ofile << "@    xaxis  ticklabel stagger 0\n";
	ofile << "@    xaxis  ticklabel place normal\n";
	ofile << "@    xaxis  ticklabel offset auto\n";
	ofile << "@    xaxis  ticklabel offset 0.000000 , 0.010000\n";
	ofile << "@    xaxis  ticklabel start type auto\n";
	ofile << "@    xaxis  ticklabel start 0.000000\n";
	ofile << "@    xaxis  ticklabel stop type auto\n";
	ofile << "@    xaxis  ticklabel stop 0.000000\n";
	ofile << "@    xaxis  ticklabel char size 2.000000\n";
	ofile << "@    xaxis  ticklabel font 0\n";
	ofile << "@    xaxis  ticklabel color 1\n";
	ofile << "@    xaxis  tick place both\n";
	ofile << "@    xaxis  tick spec type none\n";
	ofile << "@    yaxis  on\n";
	ofile << "@    yaxis  type zero false\n";
	ofile << "@    yaxis  offset 0.000000 , 0.000000\n";
	ofile << "@    yaxis  bar on\n";
	ofile << "@    yaxis  bar color 1\n";
	ofile << "@    yaxis  bar linestyle 1\n";
	ofile << "@    yaxis  bar linewidth 1.0\n";
	ofile << "@    yaxis  label \"" << yLabel << "\"\n";
	ofile << "@    yaxis  label layout para\n";
	ofile << "@    yaxis  label place auto\n";
	ofile << "@    yaxis  label char size 3.000000\n";
	ofile << "@    yaxis  label font 0\n";
	ofile << "@    yaxis  label color 1\n";
	ofile << "@    yaxis  label place normal\n";
	ofile << "@    yaxis  tick on\n";
	ofile << "@    yaxis  tick major " << TickY << "\n";
	ofile << "@    yaxis  tick minor ticks 0\n";
	ofile << "@    yaxis  tick default 6\n";
	ofile << "@    yaxis  tick place rounded true\n";
	ofile << "@    yaxis  tick in\n";
	ofile << "@    yaxis  tick major size 1.000000\n";
	ofile << "@    yaxis  tick major color 1\n";
	ofile << "@    yaxis  tick major linewidth 1.0\n";
	ofile << "@    yaxis  tick major linestyle 1\n";
	ofile << "@    yaxis  tick major grid off\n";
	ofile << "@    yaxis  tick minor color 1\n";
	ofile << "@    yaxis  tick minor linewidth 1.0\n";
	ofile << "@    yaxis  tick minor linestyle 1\n";
	ofile << "@    yaxis  tick minor grid off\n";
	ofile << "@    yaxis  tick minor size 0.500000\n";
	ofile << "@    yaxis  ticklabel on\n";
	ofile << "@    yaxis  ticklabel format general\n";
	ofile << "@    yaxis  ticklabel prec 5\n";
	ofile << "@    yaxis  ticklabel formula \"\"\n";
	ofile << "@    yaxis  ticklabel append \"\"\n";
	ofile << "@    yaxis  ticklabel prepend \"\"\n";
	ofile << "@    yaxis  ticklabel angle 0\n";
	ofile << "@    yaxis  ticklabel skip 0\n";
	ofile << "@    yaxis  ticklabel stagger 0\n";
	ofile << "@    yaxis  ticklabel place normal\n";
	ofile << "@    yaxis  ticklabel offset auto\n";
	ofile << "@    yaxis  ticklabel offset 0.000000 , 0.010000\n";
	ofile << "@    yaxis  ticklabel start type auto\n";
	ofile << "@    yaxis  ticklabel start 0.000000\n";
	ofile << "@    yaxis  ticklabel stop type auto\n";
	ofile << "@    yaxis  ticklabel stop 0.000000\n";
	ofile << "@    yaxis  ticklabel char size 2.000000\n";
	ofile << "@    yaxis  ticklabel font 0\n";
	ofile << "@    yaxis  ticklabel color 1\n";
	ofile << "@    yaxis  tick place both\n";
	ofile << "@    yaxis  tick spec type none\n";
	ofile << "@    altxaxis  off\n";
	ofile << "@    altyaxis  off\n";
	ofile << "@    legend on\n";
	ofile << "@    legend loctype view\n";
	ofile << "@    legend 0.85, 0.8\n";
	ofile << "@    legend box color 1\n";
	ofile << "@    legend box pattern 1\n";
	ofile << "@    legend box linewidth 1.0\n";
	ofile << "@    legend box linestyle 1\n";
	ofile << "@    legend box fill color 0\n";
	ofile << "@    legend box fill pattern 1\n";
	ofile << "@    legend font 0\n";
	ofile << "@    legend char size 1.000000\n";
	ofile << "@    legend color 1\n";
	ofile << "@    legend length 8\n";
	ofile << "@    legend vgap 1\n";
	ofile << "@    legend hgap 1\n";
	ofile << "@    legend invert false\n";
	ofile << "@    frame type 0\n";
	ofile << "@    frame linestyle 1\n";
	ofile << "@    frame linewidth 1.0\n";
	ofile << "@    frame color 1\n";
	ofile << "@    frame pattern 1\n";
	ofile << "@    frame background color 0\n";
	ofile << "@    frame background pattern 0\n";
	for (size_t i = 0; i<dataSets.size(); i++)
	{
		ofile << "@    s" << i << " hidden false\n";
		ofile << "@    s" << i << " type " << dataSets[i].type << "\n";
		ofile << "@    s" << i << " symbol " << dataSets[i].symbolStyle << "\n";
		ofile << "@    s" << i << " symbol size " << dataSets[i].symbolSize << "\n";
		ofile << "@    s" << i << " symbol color " << dataSets[i].symbolColor << "\n";
		ofile << "@    s" << i << " symbol pattern 1\n";
		ofile << "@    s" << i << " symbol fill color 1\n";
		ofile << "@    s" << i << " symbol fill pattern 0\n";
		ofile << "@    s" << i << " symbol linewidth 1.0\n";
		ofile << "@    s" << i << " symbol linestyle 1\n";
		ofile << "@    s" << i << " symbol char 65\n";
		ofile << "@    s" << i << " symbol char font 0\n";
		ofile << "@    s" << i << " symbol skip 0\n";
		ofile << "@    s" << i << " line type 1\n";
		ofile << "@    s" << i << " line linestyle " << dataSets[i].lineStyle << "\n";
		ofile << "@    s" << i << " line linewidth 3.0\n";
		ofile << "@    s" << i << " line color " << dataSets[i].lineColor << "\n";
		ofile << "@    s" << i << " line pattern 1\n";
		ofile << "@    s" << i << " baseline type 0\n";
		ofile << "@    s" << i << " baseline off\n";
		ofile << "@    s" << i << " dropline off\n";
		ofile << "@    s" << i << " fill type 0\n";
		ofile << "@    s" << i << " fill rule 0\n";
		ofile << "@    s" << i << " fill color 1\n";
		ofile << "@    s" << i << " fill pattern 1\n";
		ofile << "@    s" << i << " avalue off\n";
		ofile << "@    s" << i << " avalue type 2\n";
		ofile << "@    s" << i << " avalue char size 1.000000\n";
		ofile << "@    s" << i << " avalue font 0\n";
		ofile << "@    s" << i << " avalue color 1\n";
		ofile << "@    s" << i << " avalue rot 0\n";
		ofile << "@    s" << i << " avalue format general\n";
		ofile << "@    s" << i << " avalue prec 3\n";
		ofile << "@    s" << i << " avalue prepend \"\"\n";
		ofile << "@    s" << i << " avalue append \"\"\n";
		ofile << "@    s" << i << " avalue offset 0.000000 , 0.000000\n";
		ofile << "@    s" << i << " errorbar on\n";
		ofile << "@    s" << i << " errorbar place both\n";
		ofile << "@    s" << i << " errorbar color 1\n";
		ofile << "@    s" << i << " errorbar pattern 1\n";
		ofile << "@    s" << i << " errorbar size 1.000000\n";
		ofile << "@    s" << i << " errorbar linewidth 1.0\n";
		ofile << "@    s" << i << " errorbar linestyle 1\n";
		ofile << "@    s" << i << " errorbar riser linewidth 1.0\n";
		ofile << "@    s" << i << " errorbar riser linestyle 1\n";
		ofile << "@    s" << i << " errorbar riser clip off\n";
		ofile << "@    s" << i << " errorbar riser clip length 0.100000\n";
		ofile << "@    s" << i << " comment \"\"\n";
		ofile << "@    s" << i << " legend  \"" << dataSets[i].legend << "\"\n";
	}
	for (size_t i = 0; i<dataSets.size(); i++)
	{
		const std::vector<GeometryVector> & result = dataSets[i].result;
		ofile << "@target G0.S" << i << "\n";
		ofile << "@type " << dataSets[i].type << "\n";
		if (dataSets[i].type == "xy")
		{
			for (auto iter = result.begin(); iter != result.end(); iter++)
				ofile << iter->x[0] << ' ' << iter->x[1] << '\n';
		}
		else if (dataSets[i].type == "xydy")
		{
			for (auto iter = result.begin(); iter != result.end(); iter++)
				ofile << iter->x[0] << ' ' << iter->x[1] << ' ' << iter->x[2] << '\n';
		}
		else if (dataSets[i].type == "xydxdy")
		{
			for (auto iter = result.begin(); iter != result.end(); iter++)
				ofile << iter->x[0] << ' ' << iter->x[1] << ' ' << iter->x[2] << ' ' << iter->x[3] << '\n';
		}
		else
			std::cerr << "Error in gracePlot : unrecognized set type!\n";

		ofile << '&';
	}
}

std::string PlotFunction_Grace_SetType = "xy";
//presult is an array of size NumSet.
void PlotFunction_Grace(const std::vector<GeometryVector> * presult, size_t NumSet, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel, const std::vector<std::string> & legends, const std::string & Title, double MinX, double MaxX, double TickX, double MinY, double MaxY, double TickY)
{
	gracePlot plot;
	for (int i = 0; i < NumSet; i++)
		plot.addDataSet(presult[i], (legends.size()>i) ? legends[i] : "", -1, -1, -1, -1, PlotFunction_Grace_SetType, 1.0);
	plot.MaxX = MaxX;
	plot.MaxY = MaxY;
	plot.MinX = MinX;
	plot.MinY = MinY;
	plot.TickX = TickX;
	plot.TickY = TickY;
	plot.xLabel = xLabel;
	plot.yLabel = yLabel;
	plot.Title = Title;

	plot.outputFigure(OutFilePrefix);
}
void PlotFunction_Grace(const std::vector<GeometryVector> * presult, size_t NumSet, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel, const std::vector<std::string> & legends, const std::string & Title)
{
	const std::vector<GeometryVector> & result = presult[0];
	const size_t NumData=result.size();
	if(NumData==0)
	{
		std::cerr<<"Error in PlotFunction_Grace: No Data ! Exiting function \n";
		return;
	}
	double MinX = result[0].x[0], MaxX = result[0].x[0];
	double MinY = result[0].x[1], MaxY = result[0].x[1];
	for (size_t j = 0; j<NumSet; j++)
	{
		const std::vector<GeometryVector> & result = presult[j];
		const size_t NumData=result.size();
		for(size_t i=0; i<NumData; i++)
		{
			if(MinX>result[i].x[0])
				MinX=result[i].x[0];
			if(MaxX<result[i].x[0])
				MaxX=result[i].x[0];
			if(MinY>result[i].x[1])
				MinY=result[i].x[1];
			if(MaxY<result[i].x[1])
				MaxY=result[i].x[1];
		}
	}
	double TickX=std::pow(10.0, std::floor(std::log10(MaxX-MinX)-0.1761));
	if( (MaxX-MinX)/TickX<3 )
		TickX*=0.5;
	else if( (MaxX-MinX)/TickX>6 )
		TickX*=2;
 
	double TickY=std::pow(10.0, std::floor(std::log10(MaxY-MinY)-0.1761));
	if( (MaxY-MinY)/TickY<3 )
		TickY*=0.5;
	else if( (MaxY-MinY)/TickY>6 )
		TickY*=2;

	PlotFunction_Grace(presult, NumSet, OutFilePrefix, xLabel, yLabel, legends, Title, MinX, MaxX, TickX, MinY, MaxY, TickY);
}
void PlotFunction_Grace(const std::vector<GeometryVector> & result, const std::string & OutFilePrefix, const std::string & xLabel, const std::string & yLabel, const std::string & Title)
{
	PlotFunction_Grace(&result, 1, OutFilePrefix, xLabel, yLabel, std::vector<std::string>(), Title);
}

void ReadGraceData(std::vector<std::vector<GeometryVector> > & result, std::istream & ifile)
{
	if (!ifile.good())
	{
		std::cout << "Error in ReadGraceData : input file is not valid.\n";
		return;
	}
	result.clear();
	result.push_back(std::vector<GeometryVector>());
	ifile.ignore(10000, '\n');
	std::string line;
	while (!ifile.eof())
	{
		char c = ifile.peek();
		if (c == '@')
			ifile.ignore(10000, '\n');
		else if (c == '#')
			ifile.ignore(10000, '\n');
		else if (c == '&')
		{
			ifile.ignore(10000, '\n');
			result.push_back(std::vector<GeometryVector>());
		}
		else
		{
			std::getline(ifile, line);
			std::stringstream ss(line);
			GeometryVector vec;
			for (int i = 0; i<4 && ss >> vec.x[i]; i++)
				vec.SetDimension(i + 1);
			result.back().push_back(vec);
		}
	}
	result.pop_back();
}
void ReadGraceData(std::vector<std::vector<GeometryVector> > & result, const std::string & InFilePrefix)
{
	std::string name = InFilePrefix + std::string(".agr");
	std::fstream ifile(name.c_str(), std::fstream::in);
	ReadGraceData(result, ifile);
}
