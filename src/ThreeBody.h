#ifndef THREEBODY_INCLUDED
#define THREEBODY_INCLUDED

#include "PeriodicCellList.h"
#include <vector>
#include "LatticeSumSolver.h"

struct ThreeBodyThetaSeries
{
	PeriodicCellList<size_t> * pSideLengths;
	double ListSideLength;
	double Rcut;
	ThreeBodyThetaSeries(double Rc)
	{
		this->Rcut=Rc;
		this->ListSideLength=Rc*2;
		std::vector<GeometryVector> bas;
		bas.push_back(GeometryVector(this->ListSideLength, 0, 0));
		bas.push_back(GeometryVector(0, this->ListSideLength, 0));
		bas.push_back(GeometryVector(0, 0, this->ListSideLength));
		this->pSideLengths = new PeriodicCellList<size_t>(3, &bas[0], 0.1*Rc, false);
	}
	ThreeBodyThetaSeries(const ThreeBodyThetaSeries & source) : pSideLengths(new PeriodicCellList<size_t>(*source.pSideLengths)), ListSideLength(source.ListSideLength), Rcut(source.Rcut)
	{
	}
	GeometryVector ConvertToRelativeCoordinate(double * ls)//convert side lengths to RelativeCoordinate to facilitate searching in the list
	{
		std::sort(ls, ls+3);
		GeometryVector result(ls[0]/this->ListSideLength, ls[1]/this->ListSideLength, ls[2]/this->ListSideLength);
		return result;
	}
	signed long SearchTerm(double * ls, double tol)//search the nearest term with side lengths ls
		//return -1 if not found with tol distance
	{
		size_t nbr=this->pSideLengths->NumParticle();
		if(nbr==0)
			return -1;
		signed long  result=-1;
		double currentdistance2=tol*tol;
		this->pSideLengths->IterateThroughNeighbors(this->ConvertToRelativeCoordinate(ls), tol, 
			[&result, &currentdistance2](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Source)->void
		{
			if(shift.Modulus2()<currentdistance2)
			{
				result=Source;
				currentdistance2=shift.Modulus2();
			}
		});
		return result;
	}
	void AddTerm(double * ls)
	{
		this->pSideLengths->Insert(1, this->ConvertToRelativeCoordinate(ls));
	}
	~ThreeBodyThetaSeries()
	{
		delete this->pSideLengths;
	}
};


//v is the volume V, not specific volume V/N
double ThreeBodyThetaSeriesDistance(ThreeBodyThetaSeries s1, double v1, size_t n1, ThreeBodyThetaSeries s2, double v2, size_t n2, double CutOff1, DimensionType dim);
ThreeBodyThetaSeries GetThreeBodyThetaSeries(const Configuration & crystal, double Rc);
double ThreeBodySolverDistance(LatticeSumSolver * a, LatticeSumSolver * b);
#endif