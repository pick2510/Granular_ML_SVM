#include "ThreeBody.h"

//v is the volume V, not specific volume V/N
double ThreeBodyThetaSeriesDistance(ThreeBodyThetaSeries s1, double v1, size_t n1, ThreeBodyThetaSeries s2, double v2, size_t n2, double CutOff1, DimensionType dim)
{
	double Origin[3]={0.0, 0.0, 0.0};
	double scale=std::pow((v1/n1)/(v2/n2), 1.0/dim);//rescale factor
	double scale1=std::pow((v1/n1), 1.0/dim);
	double scale2=std::pow((v2/n2), 1.0/dim);
	double cutoff=s2.Rcut*scale;
	if(cutoff>s1.Rcut)
		cutoff=s1.Rcut;

	//rescale all series numbers
	for(size_t i=0; i<s1.pSideLengths->NumParticle(); i++)
	{
		//auto pa=s1.pSideLengths->GetParticle(i);
		s1.pSideLengths->GetCharacteristics(i)*=n2;
	}
	for(size_t i=0; i<s2.pSideLengths->NumParticle(); i++)
	{
		//auto pa=s2.pSideLengths->GetParticle(i);
		s2.pSideLengths->GetCharacteristics(i)*=n1;
	}

	double result=0;
	while(s2.pSideLengths->NumParticle()!=0)
	{
		signed long ShortestTerm=-1;
		for(double tol=1e-1; ;tol*=2)
		{
			ShortestTerm=s2.SearchTerm(Origin, tol);
			if(ShortestTerm!=(-1))
				break;
		}
		//auto pAtom=s2.pSideLengths->GetParticle(ShortestTerm);
		GeometryVector term=s2.pSideLengths->GetCartesianCoordinates(ShortestTerm)*scale;
		if(term.x[2]>cutoff)
		{
			s2.pSideLengths->DeleteParticle(ShortestTerm);
			continue;
		}

		//search for nearest term in s1
		signed long NearestTerm=-1;
		for(double tol=1e-3; ;tol*=2)
		{
			if(tol>cutoff)
				break;
			NearestTerm=s1.SearchTerm(term.x, tol);
			if(NearestTerm!=(-1))
				break;
		}

		if(NearestTerm==(-1))
		{
			//no nearest term remaining in s2
			break;
		}
		else
		{
			//auto pAtom2=s1.pSideLengths->GetParticle(NearestTerm);
			size_t nbr=std::min(s2.pSideLengths->GetCharacteristics(ShortestTerm), s1.pSideLengths->GetCharacteristics(NearestTerm));
			GeometryVector dist=s1.pSideLengths->GetCartesianCoordinates(NearestTerm)-s2.pSideLengths->GetCartesianCoordinates(ShortestTerm);

			double MaxDist=std::max(term.x[2], s1.pSideLengths->GetCartesianCoordinates(NearestTerm).x[2]);
			result+=nbr*dist.Modulus2()*exp((-1.0)*MaxDist/scale1)*(MaxDist-CutOff1)*(MaxDist-CutOff1);
			s2.pSideLengths->GetCharacteristics(ShortestTerm)-=nbr;
			if(s2.pSideLengths->GetCharacteristics(ShortestTerm)==0)
				s2.pSideLengths->DeleteParticle(ShortestTerm);

			s1.pSideLengths->GetCharacteristics(NearestTerm)-=nbr;
			if(s1.pSideLengths->GetCharacteristics(NearestTerm)==0)
				s1.pSideLengths->DeleteParticle(NearestTerm);
		}
	}
	result=result/scale1/scale1/n1/n2;
	return result;
}
ThreeBodyThetaSeries GetThreeBodyThetaSeries(const Configuration & crystal, double Rc)
{
	ThreeBodyThetaSeries result(Rc);

	DimensionType dim=crystal.GetDimension();
	struct neighborterm
	{
		GeometryVector shift;
		signed long periodicshift[ ::MaxDimension];
		size_t atomid;
	};
	std::vector<neighborterm> neighbors;
	for(size_t i=0; i<crystal.NumParticle(); i++)
	{
		//generate neighbor list, the list contains all neighbors that are "latter" than atom i at the original periodic cell
		//"latter" when sort atoms by:
		//1. Number
		//2. PeriodicShift
		neighbors.clear();
		crystal.IterateThroughNeighbors(i, Rc, [&neighbors, &i, &dim] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void
		{
			if(SourceAtom<i)
				return;
			if(SourceAtom==i)
			{
				bool IsLatter=false;
				for(DimensionType d=0; d<dim; d++)
				{
					if(PeriodicShift[d]>0)
					{
						IsLatter=true;
						break;
					}
					if(PeriodicShift[d]<0)
					{
						IsLatter=false;
						break;
					}
				}
				if(IsLatter==false)
					return;
			}
			//debug temp
			if(neighbors.size()==6 || neighbors.size()==14)
			{
				int ii=0;
			}

			neighborterm n;
			n.shift=shift;
			std::memcpy(n.periodicshift, PeriodicShift, ::MaxDimension*sizeof(signed long));
			n.atomid=SourceAtom;
			neighbors.push_back(n);
		});

		//generate 3-body series
		double ls[3];
		for(auto iter=neighbors.begin(); iter!=neighbors.end(); iter++)
		{
			for(auto iter2=iter+1; iter2!=neighbors.end(); iter2++)
			{
				ls[0]=std::sqrt(iter->shift.Modulus2());
				ls[1]=std::sqrt(iter2->shift.Modulus2());
				GeometryVector side3=iter->shift-iter2->shift;
				ls[2]=std::sqrt(side3.Modulus2());
				if(ls[2]>Rc)
					continue;

				signed long index=result.SearchTerm(ls, ::LengthPrecision);
				if(index==(-1))
					result.AddTerm(ls);
				else
				{
					result.pSideLengths->GetCharacteristics(index)++;
				}
			}
		}
	}
	return result;
}
const double CompareThreeBodyThetaSeriesCutoff=2.1;
double ThreeBodySolverDistance(LatticeSumSolver * a, LatticeSumSolver * b)
{
	Configuration ca=a->GetStructure();
	DimensionType dim=ca.GetDimension();
	double scalea=std::pow(ca.PeriodicVolume()/ca.NumParticle(), 1.0/dim);
	ThreeBodyThetaSeries ta=GetThreeBodyThetaSeries(ca, CompareThreeBodyThetaSeriesCutoff*scalea);
	Configuration cb=b->GetStructure();
	assert(dim==cb.GetDimension());
	double scaleb=std::pow(cb.PeriodicVolume()/cb.NumParticle(), 1.0/dim);
	ThreeBodyThetaSeries tb=GetThreeBodyThetaSeries(cb, CompareThreeBodyThetaSeriesCutoff*scaleb);

	return ThreeBodyThetaSeriesDistance(ta, ca.PeriodicVolume(), ca.NumParticle(), tb, cb.PeriodicVolume(), cb.NumParticle(), CompareThreeBodyThetaSeriesCutoff*scalea, dim);
}
