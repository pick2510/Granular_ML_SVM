#include "PhononFrequency.h"
#include <iostream>
#include <string>
#include <cassert>
#ifdef USE_MATHGL
#include <mgl2/mgl.h>
#endif

double LowestFrequencySquared;//use it to broadcast calculation results to other functions

void PlotPhononFrequencies_GenerateFigure(const std::vector<double> & results, size_t NumCurves, size_t NumX, size_t NumKPointsPerLine, std::vector<char *> & Labels, char * Prefix, double Biggest, double Smallest)
{
	#ifndef USE_MATHGL
	std::cerr<<"Error in PlotPhononFrequencies: MathGL is not enabled \n";
#else
	mglData data(&results[0], results.size()/NumCurves, NumCurves);
	mglData datax(NumX);
	for(size_t i=0; i<NumX; i++)
		datax.a[i]=static_cast<double>(i)/NumKPointsPerLine;

	//prepare ticks
	mglData datat(Labels.size());
	std::string label(" ");
	for(size_t i=0; i<Labels.size(); i++)
	{
		datat.a[i]=static_cast<double>(i);
	}

	//prepare graph
	mglGraph graph;
	if(PlotPhononFrequencies_MaxFrequencySquared_Override>0)
		Biggest=PlotPhononFrequencies_MaxFrequencySquared_Override;
	graph.SetRange('y', Smallest, Biggest);
	graph.SetRange('x', 0, NumX/NumKPointsPerLine);
	graph.SetOrigin(0, Smallest);
	graph.SetFontSize(6);
	graph.SetTicksVal('x', datat, label.c_str());
	graph.Axis();
	graph.Label('y', "\\omega^2");
			graph.SetRange('y', Smallest, Biggest);
			graph.SetRange('x', 0, NumX/NumKPointsPerLine);
			graph.SetOrigin(NumX/NumKPointsPerLine, Biggest);
			graph.SetTickTempl('x', " ");
			graph.SetTickTempl('y', " ");
			graph.Axis();
	for(size_t i=0; i<NumCurves; i++)
		graph.Plot(datax, data.SubData(i), "r");
	for(size_t i=0; i<Labels.size(); i++)
		graph.Puts(mglPoint(static_cast<mreal>(i), Biggest/(-10.0)), Labels[i]);
	for(size_t i=0; i<Labels.size(); i++)
	{
		mglData tx(2);
		tx.a[0]=i;
		tx.a[1]=i;
		mglData ty(2);
		ty.a[0]=Smallest;
		ty.a[1]=Biggest;
		graph.Plot(tx, ty, "k1-3");
	}

	std::string filename;
	filename.assign(Prefix);
	filename+="_PhononFrequencies";
	filename += FigureFormat;
	graph.WriteFrame(filename.c_str());
	if( ::AlsoWriteEPS)
	{
		filename.assign(Prefix);
		filename+="_PhononFrequencies.eps";
		graph.WriteEPS(filename.c_str());
	}

#endif
}
void PlotPhononFrequencies(const std::vector<GeometryVector> & kEndPoints, const Configuration & crystal, PairPotential & pot, char * Prefix, std::vector<char *> & Labels, double Mass, size_t NumKPointsPerLine)
{
	assert(Labels.size()==kEndPoints.size());
	if(kEndPoints.size()==0 || kEndPoints.size()==1)
	{
		std::cerr<<"Error in PlotPhononFrequencies: not enough k end points\n";
		return;
	}

	std::vector<double> results;
	size_t NumCurves;
	double Biggest=0;
	double Smallest=0;
	for(auto iter=kEndPoints.begin(); iter!=kEndPoints.end()-1; iter++)
	{
		for(size_t j=0; j<NumKPointsPerLine; j++)
		{
			std::vector<double> partresult;
			double i=static_cast<double>(j)/NumKPointsPerLine;
			GeometryVector k=(*iter)*(1.0-i)+(*(iter+1))*i;
			PhononFrequency(k, crystal, pot, partresult, Mass);
			NumCurves=partresult.size();
			for(auto iter2=partresult.begin(); iter2!=partresult.end(); iter2++)
				results.push_back(*iter2);
			if(Biggest<partresult.back())
				Biggest=partresult.back();
			if(Smallest>partresult.front())
				Smallest=partresult.front();
		}
	}
	size_t NumX=results.size()/NumCurves;
	std::vector<double> datax2(NumX);
	for(size_t i=0; i<NumX; i++)
		datax2[i]=static_cast<double>(i)/NumKPointsPerLine;

	std::string filename(Prefix);
	filename.assign(Prefix);
	filename+="_PhononFrequencies.txt";

	std::fstream ofile(filename.c_str(), std::fstream::out);
	for(auto iter=Labels.begin(); iter!=Labels.end(); iter++)
		ofile<<*iter<<" \t";
	ofile<<'\n';
	for(size_t i=0; i<datax2.size(); i++)
	{
		ofile<<datax2[i]<<" \t";
		for(size_t j=0; j<NumCurves; j++)
			ofile<<results[i*NumCurves+j]<<" \t";
		ofile<<'\n';
	}

	if(Smallest<LowestFrequencySquared)
		LowestFrequencySquared=Smallest;

	PlotPhononFrequencies_GenerateFigure(results, NumCurves, NumX, NumKPointsPerLine, Labels, Prefix, Biggest, Smallest);
}



void PlotPhononFrequencies_3DHexagonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass, size_t NumKPointsPerLine)
{
	assert(crystal.GetDimension()==3);
	std::vector<GeometryVector> ks;
	std::vector<char *> ls;
	GeometryVector gamma(0.0, 0.0, 0.0);
	GeometryVector M((crystal.GetReciprocalBasisVector(0))*0.5);
	GeometryVector K((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1))*(1.0/3));
	GeometryVector A((crystal.GetReciprocalBasisVector(2))*0.5);
	GeometryVector L=M+A;
	GeometryVector H=K+A;

	ks.push_back(gamma);
	ks.push_back(M);
	ks.push_back(K);
	ks.push_back(gamma);
	ks.push_back(A);
	ks.push_back(L);
	ks.push_back(H);
	ks.push_back(A);
	ls.push_back("\\Gamma");
	ls.push_back("M");
	ls.push_back("K");
	ls.push_back("\\Gamma");
	ls.push_back("A");
	ls.push_back("L");
	ls.push_back("H");
	ls.push_back("A");

	PlotPhononFrequencies(ks, crystal, pot, Prefix, ls, Mass, NumKPointsPerLine);
}


void PlotPhononFrequencies_3DCubic(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass, size_t NumKPointsPerLine)
{
	assert(crystal.GetDimension()==3);
	std::vector<GeometryVector> ks;
	std::vector<char *> ls;
	GeometryVector gamma(0.0, 0.0, 0.0);
	GeometryVector X((crystal.GetReciprocalBasisVector(1))*0.5);
	GeometryVector M((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1))*0.5);
	GeometryVector R((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1)+crystal.GetReciprocalBasisVector(2))*0.5);

	ks.push_back(gamma);
	ks.push_back(X);
	ks.push_back(M);
	ks.push_back(gamma);
	ks.push_back(R);
	ks.push_back(X);
	ks.push_back(M);
	ks.push_back(R);
	ls.push_back("\\Gamma");
	ls.push_back("X");
	ls.push_back("M");
	ls.push_back("\\Gamma");
	ls.push_back("R");
	ls.push_back("X");
	ls.push_back("M");
	ls.push_back("R");

	PlotPhononFrequencies(ks, crystal, pot, Prefix, ls, Mass, NumKPointsPerLine);
}
void PlotPhononFrequencies_3DTetragonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass, size_t NumKPointsPerLine)
{
	assert(crystal.GetDimension()==3);
	std::vector<GeometryVector> ks;
	std::vector<char *> ls;
	GeometryVector gamma(0.0, 0.0, 0.0);
	GeometryVector X((crystal.GetReciprocalBasisVector(1))*0.5);
	GeometryVector Z((crystal.GetReciprocalBasisVector(2))*0.5);
	GeometryVector M((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1))*0.5);
	GeometryVector A((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1)+crystal.GetReciprocalBasisVector(2))*0.5);
	GeometryVector R((crystal.GetReciprocalBasisVector(1)+crystal.GetReciprocalBasisVector(2))*0.5);

	ks.push_back(gamma);
	ks.push_back(X);
	ks.push_back(M);
	ks.push_back(gamma);
	ks.push_back(Z);
	ks.push_back(R);
	ks.push_back(A);
	ks.push_back(Z);
	ks.push_back(X);
	ks.push_back(R);
	ks.push_back(M);
	ks.push_back(A);
	ls.push_back("\\Gamma");
	ls.push_back("X");
	ls.push_back("M");
	ls.push_back("\\Gamma");
	ls.push_back("Z");
	ls.push_back("R");
	ls.push_back("A");
	ls.push_back("Z");
	ls.push_back("X");
	ls.push_back("R");
	ls.push_back("M");
	ls.push_back("A");

	PlotPhononFrequencies(ks, crystal, pot, Prefix, ls, Mass, NumKPointsPerLine);
}
void PlotPhononFrequencies_3DFCC(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass, size_t NumKPointsPerLine)
{
	assert(crystal.GetDimension()==3);
	std::vector<GeometryVector> ks;
	std::vector<char *> ls;
	GeometryVector gamma(0.0, 0.0, 0.0);
	GeometryVector X((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(2))*0.5);
	GeometryVector W((2.0*crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1)+3.0*crystal.GetReciprocalBasisVector(2))*0.25);
	GeometryVector U((5.0*crystal.GetReciprocalBasisVector(0)+4.0*crystal.GetReciprocalBasisVector(1)+5.0*crystal.GetReciprocalBasisVector(2))*0.125);
	GeometryVector L((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1)+crystal.GetReciprocalBasisVector(2))*0.5);
	GeometryVector K((3.0*crystal.GetReciprocalBasisVector(0)+3.0*crystal.GetReciprocalBasisVector(1)+6.0*crystal.GetReciprocalBasisVector(2))*0.125);

	ks.push_back(gamma);
	ks.push_back(X);
	ks.push_back(W);
	ks.push_back(K);
	ks.push_back(gamma);
	ks.push_back(L);
	ks.push_back(U);
	ks.push_back(W);
	ks.push_back(L);
	ks.push_back(K);
	ks.push_back(U);
	ks.push_back(X);
	ls.push_back("\\Gamma");
	ls.push_back("X");
	ls.push_back("W");
	ls.push_back("K");
	ls.push_back("\\Gamma");
	ls.push_back("L");
	ls.push_back("U");
	ls.push_back("W");
	ls.push_back("L");
	ls.push_back("K");
	ls.push_back("U");
	ls.push_back("X");

	PlotPhononFrequencies(ks, crystal, pot, Prefix, ls, Mass, NumKPointsPerLine);
}
void PlotPhononFrequencies_2DHexagonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass, size_t NumKPointsPerLine)
{
	assert(crystal.GetDimension()==2);
	std::vector<GeometryVector> ks;
	std::vector<char *> ls;
	GeometryVector gamma(0.0, 0.0);
	GeometryVector M((crystal.GetReciprocalBasisVector(0))*0.5);
	GeometryVector K((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1))*(1.0/3));

	ks.push_back(K);
	ks.push_back(gamma);
	ks.push_back(M);
	ls.push_back("K");
	ls.push_back("\\Gamma");
	ls.push_back("M");

	PlotPhononFrequencies(ks, crystal, pot, Prefix, ls, Mass, NumKPointsPerLine);
}
void PlotPhononFrequencies_2DCubic(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass, size_t NumKPointsPerLine)
{
	assert(crystal.GetDimension()==2);
	std::vector<GeometryVector> ks;
	std::vector<char *> ls;
	GeometryVector gamma(0.0, 0.0);
	GeometryVector X((crystal.GetReciprocalBasisVector(1))*0.5);
	GeometryVector M((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1))*0.5);

	ks.push_back(gamma);
	ks.push_back(M);
	ks.push_back(X);
	ks.push_back(gamma);
	ls.push_back("\\Gamma");
	ls.push_back("M");
	ls.push_back("X");
	ls.push_back("\\Gamma");

	PlotPhononFrequencies(ks, crystal, pot, Prefix, ls, Mass, NumKPointsPerLine);
}
void PlotPhononFrequencies_2DRectangle(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass, size_t NumKPointsPerLine)
{
	assert(crystal.GetDimension()==2);
	std::vector<GeometryVector> ks;
	std::vector<char *> ls;
	GeometryVector gamma(0.0, 0.0);
	GeometryVector X((crystal.GetReciprocalBasisVector(0))*0.5);
	GeometryVector Y((crystal.GetReciprocalBasisVector(1))*0.5);
	GeometryVector S((crystal.GetReciprocalBasisVector(0)+crystal.GetReciprocalBasisVector(1))*0.5);

	ks.push_back(gamma);
	ks.push_back(X);
	ks.push_back(S);
	ks.push_back(Y);
	ks.push_back(gamma);
	ls.push_back("\\Gamma");
	ls.push_back("X");
	ls.push_back("S");
	ls.push_back("Y");
	ls.push_back("\\Gamma");

	PlotPhononFrequencies(ks, crystal, pot, Prefix, ls, Mass, NumKPointsPerLine);
}