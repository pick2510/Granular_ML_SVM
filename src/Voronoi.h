#ifndef VORONOI_INCLUDED
#define VORONOI_INCLUDED

#include <set>
#include "PeriodicCellList.h"

//calculate all vertexes of the voronoi cell of particle i, in vertices
//if pLinkage is not nullptr, it will be filled with linkage information of the format:
// (*pLinkage)[i*(Dimension+1)+j] is the index of a vertex linked with vertex i
// (*pLinkage)[i*(Dimension+1)+j] is -1 if vertex i have less than j linked vertexes.
void GetVoronoiCell(const PeriodicCellList<Empty> & c, size_t i, std::vector<GeometryVector> & vertices, std::vector<signed long> * pLinkage = nullptr, std::set<size_t> * pNeighbors = nullptr);

double PolygonVolume(const std::vector<GeometryVector> & Vertices, DimensionType d);
void PlotVoronoi2D(std::string prefix, const Configuration & List, const std::string & remark="");

//Num=-1 : display voronoi of all particles Num>=0 : display voronoi of a specific particle
void DisplayConfiguration_3D_WithVoronoi(const Configuration & list, double Radius=0.0, long Num=-1);


//Probability density function of the voronoi cell volume
//clear Result, then fill it with the PDF calculated from Config
//Result is filled with 4-dimensional GeometryVectors, the elements are:
//(r, p(r), \delta r, \delta p(r) )
//gets SampleSize number of sample pair distances, use it to generate bins, then count all pair distances
//ResolutionPreference : >1 to get better r resolution, <1 to get better g_2 resolution
void VoronoiVolumeDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleSize=500000, double ResolutionPreference=1.0);


void VoronoiNumSidesDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig, std::vector<GeometryVector> & Result);

void VoronoiVolumeCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, double MaxDistance, double RPrecision);


#endif