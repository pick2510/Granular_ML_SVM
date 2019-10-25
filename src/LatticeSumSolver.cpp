#include "LatticeSumSolver.h"
#include "etc.h"

#include <algorithm>
#include <vector>


void GetLatticeTerm(Configuration & crystal, std::vector<LatticeTerm> & result, double Rc)
{
	result.clear();

	std::vector<double> distance2s;
	for(size_t i=0; i<crystal.NumParticle(); i++)
	{
		//Configuration::particle * temp = crystal.GetParticle(i);
		crystal.IterateThroughNeighbors(i, Rc, [&distance2s] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) -> void 
		{
			distance2s.push_back(shift.Modulus2());
		});
	}

	if(distance2s.size()>crystal.NumParticle())
	{
		std::stable_sort(distance2s.begin(), distance2s.end());

		double nowr2=0;
		unsigned long count=0;
		for(auto iterd=distance2s.begin(); iterd!=distance2s.end(); iterd++)
		{
			if(std::fabs(nowr2-*iterd)<LengthPrecision)
				count++;
			else
			{
				LatticeTerm temp;
				temp.number=count;
				temp.distance=std::sqrt(nowr2);
				result.push_back(temp);

				nowr2=*iterd;
				count=1;
			}
		}
		LatticeTerm temp;
		temp.number=count;
		temp.distance=std::sqrt(nowr2);
		result.push_back(temp);
	}
}


LatticeSumSolver::LatticeSumSolver()
{
	this->CurrentRc=0;
	this->OriginalDensity=1;
	this->Dimension=1;
	this->NumParticle=0;
}
LatticeSumSolver::~LatticeSumSolver()
{
}
void LatticeSumSolver::UpdateTerms(double NewRc)
{
	Configuration temp=this->GetStructure();
	this->Terms.clear();
	GetLatticeTerm(temp, Terms, NewRc);
	this->NumParticle=temp.NumParticle();
	this->CurrentRc=NewRc;
	this->Dimension=temp.GetDimension();
	this->OriginalDensity=temp.NumParticle()/temp.PeriodicVolume();
}
const double RescalculateFactor=1.1;
double LatticeSumSolver::LatticeSum(double Volume, IsotropicPairPotential & potential)
{
	//rescale Rc using density
	Configuration temp=this->GetStructure();
	this->OriginalDensity=temp.NumParticle()/temp.PeriodicVolume();
	double a=std::pow(this->OriginalDensity*Volume, 1.0/this->Dimension);//rescale factor
	double Rc=std::sqrt(potential.Rcut2)/a;
	if(this->CurrentRc<Rc)
	{
		this->UpdateTerms(RescalculateFactor*Rc);
		return this->LatticeSum(Volume, potential);
	}

	double result=0;

	double tempoffset;
	if(this->Terms.size()!=0)
	{
		for(std::vector<LatticeTerm>::iterator iter=this->Terms.begin()+1; iter!=this->Terms.end(); iter++)
		{
			if(iter->distance>Rc)break;//both iter->distance and Rc are re-scaled
			tempoffset=iter->distance*a;
			result+=iter->number*(potential.IsotropicPairEnergy(tempoffset, "", ""));
		}
	}

	result /= this->NumParticle;

	result /= 2;

	return result;
}

bool SameSolver(LatticeSumSolver * a, LatticeSumSolver * b)
{
	return (SolverDistance(a, b)<1e-7);
}

double SolverDistance(LatticeSumSolver * a, LatticeSumSolver * b)
{
	Configuration l1=a->GetStructure();
	Configuration l2=b->GetStructure();
	DimensionType dim=l2.GetDimension();
	if(l1.GetDimension()!=dim)
		return false;
	double ascale=std::pow(l1.PeriodicVolume()/l1.NumParticle(), 1.0/dim);//rescale factor
	double bscale=std::pow(l2.PeriodicVolume()/l2.NumParticle(), 1.0/dim);//rescale factor

	if(a->CurrentRc<TwoBodyDistance_MaxLength*ascale)
		a->UpdateTerms(TwoBodyDistance_MaxLength*ascale);
	if(b->CurrentRc<TwoBodyDistance_MaxLength*bscale)
		b->UpdateTerms(TwoBodyDistance_MaxLength*bscale);

	//rescale everything
	std::vector<LatticeTerm> RescaledA, RescaledB;
	for(auto iter=a->Terms.begin(); iter!=a->Terms.end(); iter++)
	{
		LatticeTerm temp;
		temp.distance=iter->distance/ascale;
		temp.number=iter->number*l2.NumParticle();
		RescaledA.push_back(temp);
	}
	for(auto iter=b->Terms.begin(); iter!=b->Terms.end(); iter++)
	{
		LatticeTerm temp;
		temp.distance=iter->distance/bscale;
		temp.number=iter->number*l1.NumParticle();
		RescaledB.push_back(temp);
	}

	std::vector<LatticeTerm>::iterator itera=RescaledA.begin(), iterb=RescaledB.begin();
	double result=0;
	while(itera!=RescaledA.end() && iterb!=RescaledB.end())
	{
		auto NumComparing=std::min(itera->number, iterb->number);

		double dist=itera->distance-iterb->distance;
		result+=NumComparing*dist*dist*std::exp((-1)*itera->distance);
		itera->number-=NumComparing;
		iterb->number-=NumComparing;

		if(itera->number==0)
			itera++;
		if(iterb->number==0)
			iterb++;
	}
	result /=(l1.NumParticle()*l2.NumParticle());

	return result;
}