#ifndef PHONONFREQUENCY_INCLUDED
#define PHONONFREQUENCY_INCLUDED

#include "GeometryVector.h"
#include "Potential.h"
#include "PeriodicCellList.h"
#include <vector>
void PhononFrequency(const GeometryVector & k, const Configuration & crystal, PairPotential & pot, std::vector<double> & result2, double Mass=1.0);
void PlotPhononFrequencies(const std::vector<GeometryVector> & kEndPoints, const Configuration & crystal, PairPotential & pot, char * Prefix, std::vector<char *> & Labels, double Mass=1.0, size_t NumKPointsPerLine=30);


extern double LowestFrequencySquared;//use it to broadcast calculation results to other functions

void PlotPhononFrequencies_3DHexagonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=30);
void PlotPhononFrequencies_3DTetragonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=30);
void PlotPhononFrequencies_3DCubic(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=30);
void PlotPhononFrequencies_2DHexagonal(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=90);
void PlotPhononFrequencies_2DCubic(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=90);
void PlotPhononFrequencies_2DRectangle(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=60);
void PlotPhononFrequencies_3DFCC(const Configuration & crystal, PairPotential & pot, char * Prefix, double Mass=1.0, size_t NumKPointsPerLine=20);

#endif