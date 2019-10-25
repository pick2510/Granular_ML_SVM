#ifndef POLYNOMIAL_INCLUDED
#define POLYNOMIAL_INCLUDED

#include "etc.h"
#include "Potential.h"
#include "ParameterSet.h"



const double Polynomial_CoreSize_Over_NeighborDistance = 0.8;
class GaussianCorePotential : public IsotropicPairPotential
{
public:
	double * data;
	size_t Order;
	double NeighborDistance;
	double a0;//0th order polynomial term
	GaussianCorePotential(DimensionType dimension, size_t Order, double RCut, double NeighborDistance, double a0) : IsotropicPairPotential(dimension, RCut), NeighborDistance(NeighborDistance), a0(a0)
	{
		this->Order=Order;
		this->data = new double [Order];
	}
	~GaussianCorePotential()
	{
		delete [] this->data;
	}
	virtual double IsotropicEnergyFormula(double distance, const AtomInfo & a1, const AtomInfo & a2)
	{
		double & RC = this->data[0];
		if(distance>RC) return 0;

		double result = 0;
		result += std::pow((this->NeighborDistance*Polynomial_CoreSize_Over_NeighborDistance/distance), 12);//hard core
		double PolyValue=this->data[Order-1];
		for(size_t i=Order-2; i>=2; i--)
		{
			PolyValue=PolyValue*distance+this->data[i];
		}
		PolyValue=PolyValue*distance;
		//a0 is fixed so that polynomial at 1st neighbor distance is 1
		PolyValue+=this->a0;

		result+=PolyValue;
		result *= (RC - distance);
		result *= (RC - distance);
		result *= std::exp((-1)*this->data[1]*distance*distance);
		return result;
	}
	virtual double IsotropicForceFormula(double distance, const AtomInfo & a1, const AtomInfo & a2)
	{
		double & RC = this->data[0];
		if(distance>RC) return 0;

		double part1 = 0;
		part1 += std::pow((this->NeighborDistance*Polynomial_CoreSize_Over_NeighborDistance/distance), 12);//hard core
		double PolyValue=this->data[Order-1];
		double dPolyValue=(Order-2)*this->data[Order-1];
		for(size_t i=Order-2; i>=2; i--)
		{
			PolyValue=PolyValue*distance+this->data[i];
			dPolyValue=dPolyValue*distance+this->data[i]*(i-1);
		}
		PolyValue=PolyValue*distance;
		//a0 is fixed so that polynomial at 1st neighbor distance is 1
		PolyValue+=this->a0;

		part1+=PolyValue;


		double part2=(RC - distance)*(RC - distance);


		double part3= std::exp((-1)*this->data[1]*distance*distance);

		double dpart1=0;
		dpart1 += std::pow((this->NeighborDistance*Polynomial_CoreSize_Over_NeighborDistance/distance), 12)*(-12.0)/distance;
		dpart1+=dPolyValue;
		double dpart2=2.0*(distance-RC);
		double dpart3=part3*(-2.0)*distance*this->data[1];

		double result=(-1.0)*(part1*part2*dpart3+part1*dpart2*part3+dpart1*part2*part3);

		//assert(std::abs((this->IsotropicPairPotential::IsotropicForceFormula(distance)-result)/result)<1e-6);
		return result;
	}
};

const double Polynomial_Max_Rcut_Over_NeighborDistance = 6.4;
const double Polynomial_MinDistance_Over_NeighborDistance = 0.3;
class PolynomialParameterization : public ParameterSet
{
public:
	size_t Order;
	DimensionType dimension;
	double NeighborDistance;
	PolynomialParameterization(DimensionType dimension, double NeighborDistance) : dimension(dimension)
	{
		this->Order=2;
		this->NeighborDistance=NeighborDistance;
		this->Parameters = new double [Order];
		for(size_t i=2; i<Order; i++)
			this->Parameters[i] =2;
		this->Parameters[0]=1.5*NeighborDistance;
		this->Parameters[1]=0.0;
	}
	PolynomialParameterization(LatticeSumSolver * Target) : dimension(Target->GetStructure().GetDimension())
	{
		this->Order=2;

		//make sure Target has some terms
		double tRc=1;
		while(Target->Terms.size()<2)
		{
			Target->UpdateTerms(tRc);
			tRc*=2;
		}
		this->NeighborDistance = Target->Terms[1].distance;

		this->Parameters = new double [Order];
		for(size_t i=2; i<Order; i++)
			this->Parameters[i] =2;
		this->Parameters[0]=Target->Terms[1].distance*1.1;
		this->Parameters[1]=0.0;
	}
	PolynomialParameterization(const PolynomialParameterization & src) : ParameterSet(src), dimension(src.dimension), Order(src.Order), NeighborDistance(src.NeighborDistance)
	{
	}
	virtual double MaxDistance(void) const
	{
		return this->NeighborDistance*Polynomial_Max_Rcut_Over_NeighborDistance;
	}
	virtual double MinDistance(void) const
	{
		return this->NeighborDistance*Polynomial_MinDistance_Over_NeighborDistance;
	}
	virtual IsotropicPairPotential * GetPotential(void) const
	{
		double PolyValueAtNeighbor=this->Parameters[Order-1];
		for(size_t i=Order-2; i>=2; i--)
		{
			PolyValueAtNeighbor=PolyValueAtNeighbor*this->NeighborDistance+this->Parameters[i];
		}
		PolyValueAtNeighbor=PolyValueAtNeighbor*this->NeighborDistance;
		//a0 is fixed so that polynomial at 1st neighbor distance is 1
		double a0=1.0-PolyValueAtNeighbor;
		GaussianCorePotential * result = new GaussianCorePotential(this->dimension, this->Order, this->Parameters[0], this->NeighborDistance, a0);
		std::memcpy(result->data, this->Parameters, sizeof(double)*Order);
		//if( ::UseCorrectedPotential)
		//	return new CorrectedIsotropicPairPotential(result);
		//else
			return result;
	}
	virtual double Penalty(void) const
	{
		return std::pow(this->Parameters[0], dimension);
		//return std::exp(15.0*this->Parameters[0]);
	}
	virtual void GetBound(double * LowerBound, double * UpperBound) const
	{
		for(size_t i=2; i<Order; i++)
		{
			LowerBound[i]=-20;
			UpperBound[i]=20;
		}
		LowerBound[0]=this->NeighborDistance*1.05;
		UpperBound[0]=this->NeighborDistance*Polynomial_Max_Rcut_Over_NeighborDistance;
		LowerBound[1]=0.0;
		UpperBound[1]=0.0;
	}
	virtual bool Evolve(void)
	{
		if(this->Order>100)
			return false;
		double * newPara = new double [this->Order+1];
		for(size_t i=0; i<this->Order; i++)
			newPara[i]=this->Parameters[i];
		newPara[this->Order]=2.0;
		delete [] this->Parameters;
		this->Parameters=newPara;
		this->Order++;

		return true;
	}
	virtual bool Iterate(void)
	{
		return false;
	}
	virtual ParameterSet * clone() const
	{
		return new PolynomialParameterization(*this);
	}
	virtual void Write(std::ostream & out) const
	{
		out<<this->dimension<<" \t";
		out<<this->Order<<" \t";
		out<<this->NeighborDistance<<" \t";
		for(size_t i=0; i<Order; i++)
			out<<this->Parameters[i]<<" \t";
		out<<'\n';

		//output a intuitive expression
		double PolyValueAtNeighbor=this->Parameters[Order-1];
		for(size_t i=Order-2; i>=2; i--)
		{
			PolyValueAtNeighbor=PolyValueAtNeighbor*this->NeighborDistance+this->Parameters[i];
		}
		PolyValueAtNeighbor=PolyValueAtNeighbor*this->NeighborDistance;
		//a0 is fixed so that polynomial at 1st neighbor distance is 1
		double a0=1.0-PolyValueAtNeighbor;
		double coresigma=this->NeighborDistance*Polynomial_CoreSize_Over_NeighborDistance;
		out<<"(("<<coresigma<<"/r)^12";
		if(a0>=0)
			out<<'+';
		out<<a0;
		for(size_t i=2; i<Order; i++)
		{
			if(this->Parameters[i]>0)
				out<<'+';
			out<<this->Parameters[i]<<"*r^"<<i-1;
		}
		out<<")*("<<this->Parameters[0]<<"-r)^2*exp(-"<<this->Parameters[1]<<"*r^2)\n";
	}
	virtual void Read(std::istream & in)
	{
		in>>this->dimension;
		in>>this->Order;
		in>>this->NeighborDistance;
		delete [] this->Parameters;
		this->Parameters = new double [this->Order];
		for(size_t i=0; i<Order; i++)
			in>>this->Parameters[i];
		std::string temp;
		std::getline(in, temp);
	}
	virtual double RecommendedTolerence(void) const
	{
		return 1e-6;
	}
	virtual double RecommendedMaxEnergy(void) const
	{
		return 1e8;
	}
	virtual size_t NumParameters(void) const
	{
		return Order;
	}
	virtual double SharpestRelativeSize(void)
	{
		return 0.2;
	}
};

const double GaussianCore_Max_Exponent_At_First_Neighbor = 3;
class GaussianCoreParameterization : public PolynomialParameterization
{
public:
	GaussianCoreParameterization(DimensionType dimension, double MinDistance) : PolynomialParameterization(dimension, MinDistance/Polynomial_MinDistance_Over_NeighborDistance)
	{
	}
	GaussianCoreParameterization(LatticeSumSolver * Target) : PolynomialParameterization(Target)
	{
	}
	GaussianCoreParameterization(const GaussianCoreParameterization & src) : PolynomialParameterization(src)
	{
	}
	virtual void GetBound(double * LowerBound, double * UpperBound) const
	{
		for(size_t i=2; i<Order; i++)
		{
			LowerBound[i]=-20;
			UpperBound[i]=20;
		}
		LowerBound[0]=this->NeighborDistance*1.05;
		UpperBound[0]=this->NeighborDistance*Polynomial_Max_Rcut_Over_NeighborDistance;
		LowerBound[1]=0.0;
		UpperBound[1]=GaussianCore_Max_Exponent_At_First_Neighbor/this->NeighborDistance/this->NeighborDistance;
	}
	virtual ParameterSet * clone() const
	{
		return new GaussianCoreParameterization(*this);
	}
};

class SoftCoreParameterization : public GaussianCoreParameterization
{
public:
	SoftCoreParameterization(DimensionType dimension, double MinDistance) : GaussianCoreParameterization(dimension, MinDistance)
	{
	}
	SoftCoreParameterization(LatticeSumSolver * Target) : GaussianCoreParameterization(Target)
	{
	}
	SoftCoreParameterization(const GaussianCoreParameterization & src) : GaussianCoreParameterization(src)
	{
	}
	virtual IsotropicPairPotential * GetPotential(void) const
	{
		double PolyValueAtNeighbor=this->Parameters[Order-1];
		for(size_t i=Order-2; i>=2; i--)
		{
			PolyValueAtNeighbor=PolyValueAtNeighbor*this->NeighborDistance+this->Parameters[i];
		}
		PolyValueAtNeighbor=PolyValueAtNeighbor*this->NeighborDistance;
		//a0 is fixed so that polynomial at 1st neighbor distance is 1
		double a0=1.0-PolyValueAtNeighbor;
		//parameter 4=0.0 will disable the hard core
		GaussianCorePotential * result = new GaussianCorePotential(this->dimension, this->Order, this->Parameters[0], 0.0, a0);
		std::memcpy(result->data, this->Parameters, sizeof(double)*Order);
		return result;
	}
	virtual void Write(std::ostream & out) const
	{
		out<<this->dimension<<" \t";
		out<<this->Order<<" \t";
		out<<this->NeighborDistance<<" \t";
		for(size_t i=0; i<Order; i++)
			out<<this->Parameters[i]<<" \t";
		out<<'\n';

		//output a intuitive expression
		double PolyValueAtNeighbor=this->Parameters[Order-1];
		for(size_t i=Order-2; i>=2; i--)
		{
			PolyValueAtNeighbor=PolyValueAtNeighbor*this->NeighborDistance+this->Parameters[i];
		}
		PolyValueAtNeighbor=PolyValueAtNeighbor*this->NeighborDistance;
		//a0 is fixed so that polynomial at 1st neighbor distance is 1
		double a0=1.0-PolyValueAtNeighbor;
		double coresigma=0;
		out<<"(("<<coresigma<<"/r)^12";
		if(a0>0)
			out<<'+';
		out<<a0;
		for(size_t i=2; i<Order; i++)
		{
			if(this->Parameters[i]>0)
				out<<'+';
			out<<this->Parameters[i]<<"*r^"<<i-1;
		}
		out<<")*("<<this->Parameters[0]<<"-r)^2*exp(-"<<this->Parameters[1]<<"*r^2)\n";
	}
};

#endif
