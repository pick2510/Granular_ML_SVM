#ifndef WITHGAUSSIANINV_INCLUDED
#define WITHGAUSSIANINV_INCLUDED


#include "ParameterSet.h"
#include "LatticeSumSolver.h"
#include "Potential.h"
#include <vector>
#include <cmath>

struct InvTerm
{
	double intensity;
	double power;
	InvTerm(double Intensity, double Power)
	{
		this->intensity=Intensity;
		this->power=(-1)*std::abs(Power);
	}
};

struct ExpTerm
{
	double intensity, alpha;
	ExpTerm(double Intensity, double Alpha)
	{
		this->intensity=Intensity;
		this->alpha=(-1)*std::abs(Alpha);
	}
};

struct GaussianWell
{
	double location, width, intensity;
	GaussianWell(double Intensity, double Width, double Location):intensity(Intensity), width(Width), location(Location)
	{
	}
};

class withGaussianInvPotential: public IsotropicPairPotential
{
public:
	std::vector<InvTerm> invterms;
	std::vector<GaussianWell> wells;
	std::vector<ExpTerm> expterms;
	double Ecut;
	withGaussianInvPotential(DimensionType Dimension, double Rcut):IsotropicPairPotential(Dimension, Rcut)
	{
		this->Ecut=0;
		this->Ecut=this->IsotropicEnergyFormula(Rcut, "", "");
	}
	virtual double IsotropicEnergyFormula(double distance, const AtomInfo & a1, const AtomInfo & a2)
	{
		double result=0;
		for(auto iter=invterms.begin(); iter!=invterms.end(); iter++)
			result+=std::pow(distance, iter->power)*iter->intensity;
		/*
		for(auto iter=this->expterms.begin(); iter!=this->expterms.end(); iter++)
		result+=iter->intensity*std::exp(iter->alpha*distance);
		*/
		for(auto iter=wells.begin(); iter!=wells.end(); iter++)
		{
			double temp=(distance-iter->location)/iter->width;
			result+=iter->intensity*std::exp((-1)*temp*temp);
		}
		return result-this->Ecut;
	}
};


const double Umin=-1;


const unsigned short MaxStage=7;
class withGaussianInvParameters: public ParameterSet
{
	int exp1;
	int exp2;
	signed short stage;
	std::vector<double> wells;
	double RCut;
	::DimensionType dimension;
	bool Global;
public:
	withGaussianInvParameters(LatticeSumSolver * TargetStructure)
	{
		this->exp1=12;
		this->exp2=6;
		this->stage=0;
		this->RCut=3.0;
		this->Parameters=new double [3];

		this->Parameters[0]=1.0;
		this->Parameters[1]=3.0;
		this->Parameters[2]=3.0;

		this->dimension=TargetStructure->GetStructure().GetDimension();
		while(TargetStructure->Terms.size()< ::MaxStage+1)
			TargetStructure->UpdateTerms(TargetStructure->CurrentRc*2+1);
		for(int i=0; i< ::MaxStage; i++)
			this->wells.push_back(TargetStructure->Terms[i+1].distance);
	}
	withGaussianInvParameters(const withGaussianInvParameters & source):ParameterSet(source), wells(source.wells)
	{
		this->exp1=source.exp1;
		this->exp2=source.exp2;
		this->stage=source.stage;
		this->dimension=source.dimension;
		this->RCut=source.RCut;
	}
	withGaussianInvParameters(withGaussianInvParameters && source):ParameterSet(source), wells(source.wells)
	{
		this->exp1=source.exp1;
		this->exp2=source.exp2;
		this->stage=source.stage;
		this->dimension=source.dimension;
		this->RCut=source.RCut;
	}

	virtual IsotropicPairPotential * GetPotential(void) const
	{
		withGaussianInvPotential * result = new withGaussianInvPotential(this->dimension, RCut);

		double rmin=this->Parameters[0];
		double dexp1=static_cast<double>(this->exp1);
		double dexp2=static_cast<double>(this->exp2);
		double k=std::pow(dexp2/dexp1, static_cast<double>(1)/(dexp1-dexp2));
		double b=Umin*std::pow(k*rmin, dexp2)/(std::pow(k, dexp1)-std::pow(k, dexp2));
		double a=std::pow(k*rmin, dexp1-dexp2)*b;

		result->invterms.push_back(InvTerm(a, dexp1));
		result->invterms.push_back(InvTerm((-1)*b, dexp2));

		result->expterms.push_back(ExpTerm(this->Parameters[1], this->Parameters[2]));

		for(size_t thiswell=0; thiswell<this->stage; thiswell++)
			result->wells.push_back(GaussianWell(this->Parameters[3*thiswell+3], 1/std::exp(this->Parameters[3*thiswell+4]), this->Parameters[3*thiswell+5]));

		//if( ::UseCorrectedPotential)
		//	return new CorrectedIsotropicPairPotential(result);
		//else
			return result;
	}
	virtual double Penalty(void) const//calculate the penalty from OptimizableParameters
	{
		if(this->Global)//disable penalty for global optimizations
			return 1;

		double result=0;

		for(size_t thiswell=0; thiswell<this->stage; thiswell++)
		{
			result+=std::abs(this->Parameters[3*thiswell+3])*this->Parameters[3*thiswell+5]*this->Parameters[3*thiswell+5]*10;
			result+=std::abs(this->Parameters[3*thiswell+4]);
		}
		result/=5;

		return result;
	}
	virtual void GetBound(double * LowerBound, double * UpperBound) const//get range of each parameters
	{
		LowerBound[0]=0.4;
		UpperBound[0]=1.5;
		LowerBound[1]=0;
		UpperBound[1]=9;
		LowerBound[2]=0.8;
		UpperBound[2]=9;
		for(int i=0; i<this->stage; i++)
		{
			LowerBound[3*i+3]=-2.0;
			UpperBound[3*i+3]=2.0;
			LowerBound[3*i+4]=-1.0;
			UpperBound[3*i+4]=3.0;
			LowerBound[3*i+5]=0.4;
			UpperBound[3*i+5]=0.8*this->RCut;
		}
	}
	virtual bool Evolve(void)//evolve to a more complex Parameterization, return true if succed, return false if can't evolve any more
	{
		if(this->stage==::MaxStage)
			return false;
		this->stage++;
		double * newParameters = new double [3*(this->stage+1)];
		std::memcpy(newParameters, this->Parameters, sizeof(double)*3*this->stage);
		newParameters[3*this->stage]=-1.0;
		newParameters[3*this->stage+1]=2.3;
		newParameters[3*this->stage+2]=this->wells[this->stage-1];
		if(newParameters[3*this->stage+2]>1.9)
			newParameters[3*this->stage+2]=1.9;

		delete [] this->Parameters;
		this->Parameters=newParameters;
		return true;
	}
	virtual bool Iterate(void)//for some non-continuous parameter, iterate through it. Intended to be used when it has never evolved
	{
		if(this->exp1==12 && this->exp2==10)
			return false;
		this->exp2+=2;
		if(this->exp2>=this->exp1)
		{
			this->exp2=4;
			this->exp1+=2;
		}
		return true;
	}
	virtual ParameterSet * clone() const
	{
		return new withGaussianInvParameters(*this);
	}
	virtual void Write(std::ostream & out) const
	{
		out<<this->exp1<<" \t"<<this->exp2<<'\n';
		out<<this->stage<<" \t"<<this->RCut<<'\n';
		out<<this->Parameters[0]<<'\n';
		out<<this->Parameters[1]<<" \t"<<this->Parameters[2]<<'\n';
		int end=3*(this->stage+1);
		for(int i=3; i<end; i++)
			out<<this->Parameters[i]<<" \t";
		out<<'\n';
	}
	virtual void Read(std::istream & in)
	{
		in>>this->exp1;
		in>>this->exp2;
		in>>this->stage;
		in>>this->RCut;
		int end=3*(this->stage+1);
		delete [] this->Parameters;
		this->Parameters = new double [end];
		for(int i=0; i<end; i++)
			in>>this->Parameters[i];
	}
	virtual double MinDistance(void) const//when r<MinDistance, the potential should have no interesting behavior, guaranteed to fix unchanged during the life of this object
	{
		return 0.3;
	}
	virtual double MaxDistance(void) const//the Rcut of the potential should not exceed this, guaranteed to fix unchanged during the life of this object
	{
		return this->RCut;
	}
	virtual size_t NumParameters(void) const
	{
		return 3*(this->stage+1);
	}
	virtual void SignalGlobalOptimization(void)//a signal that global optimization will be carried out
	{
		this->Global=true;
	}
	virtual void SignalLocalOptimization(void)//vice versa
	{
		this->Global=false;
	}
	virtual double SharpestRelativeSize(void)
	{
		return 0.02;
	}
};

const size_t TunnelNumWells=6;
const double RCut=2.5;
const double DefaultWellLocations[TunnelNumWells]={1.0000, 1.34164, 1.48324, 1.67332, 1.89737, 2};
class TunnelwithGaussianInvParameters: public ParameterSet
{
	int exp1;
	int exp2;
	signed short stage;//0: optimize for intensity, 1: optimize also for well width and location
public:
	TunnelwithGaussianInvParameters()
	{
		this->exp1=12;
		this->exp2=6;
		this->stage=0;
		this->Parameters=new double [3+::TunnelNumWells];

		this->Parameters[0]=1.0;
		this->Parameters[1]=0.0;
		this->Parameters[2]=3.0;

		for(int i=3; i<3+::TunnelNumWells; i++)
		{
			this->Parameters[i]=-1.0;
		}
	}
	TunnelwithGaussianInvParameters(const TunnelwithGaussianInvParameters & source)
	{
		this->exp1=source.exp1;
		this->exp2=source.exp2;
		this->stage=source.stage;
	}
	TunnelwithGaussianInvParameters(TunnelwithGaussianInvParameters && source)
	{
		this->exp1=source.exp1;
		this->exp2=source.exp2;
		this->stage=source.stage;
	}
	virtual IsotropicPairPotential * GetPotential(void) const
	{
		withGaussianInvPotential * result = new withGaussianInvPotential(3, RCut);

		double rmin=this->Parameters[0];
		double dexp1=static_cast<double>(this->exp1);
		double dexp2=static_cast<double>(this->exp2);
		double k=std::pow(dexp2/dexp1, static_cast<double>(1)/(dexp1-dexp2));
		double b=Umin*std::pow(k*rmin, dexp2)/(std::pow(k, dexp1)-std::pow(k, dexp2));
		double a=std::pow(k*rmin, dexp1-dexp2)*b;

		result->invterms.push_back(InvTerm(a, dexp1));
		result->invterms.push_back(InvTerm((-1)*b, dexp2));

		result->expterms.push_back(ExpTerm(this->Parameters[1], this->Parameters[2]));

		if(this->stage==0)
		{
			for(size_t thiswell=0; thiswell< ::TunnelNumWells; thiswell++)
				result->wells.push_back(GaussianWell(this->Parameters[thiswell+3], 0.1, DefaultWellLocations[thiswell]));
		}
		else if (this->stage==1)
		{
			for(size_t thiswell=0; thiswell< ::TunnelNumWells; thiswell++)
				result->wells.push_back(GaussianWell(this->Parameters[3*thiswell+3], 1/std::exp(this->Parameters[3*thiswell+4]), this->Parameters[3*thiswell+5]));
		}
		else
		{
			throw "error in TunnelwithGaussianInvParameters::GetPotential : unexpected this->stage";
		}

		return result;
	}
	virtual double Penalty(void) const//calculate the penalty from OptimizableParameters
	{
		return 0;
	}
	virtual void GetBound(double * LowerBound, double * UpperBound) const//get range of each parameters
	{
		LowerBound[0]=0.4;
		UpperBound[0]=1.0;
		LowerBound[1]=0;
		UpperBound[1]=0;
		LowerBound[2]=3;
		UpperBound[2]=3;
		if(this->stage==0)
		{
			for(int i=0; i< ::TunnelNumWells; i++)
			{
				LowerBound[i+3]=-2.0;
				UpperBound[i+3]=1.0;
			}
		}
		else
		{
			for(int i=0; i< ::TunnelNumWells; i++)
			{
				LowerBound[3*i+3]=-2.0;
				UpperBound[3*i+3]=0.0;
				LowerBound[3*i+4]=-1.0;
				UpperBound[3*i+4]=3.0;
				LowerBound[3*i+5]=::DefaultWellLocations[i]*0.9;
				UpperBound[3*i+5]=::DefaultWellLocations[i]*1.1;
			}
		}
	}
	virtual bool Evolve(void)//evolve to a more complex parametrization, return true if succed, return false if can't evolve any more
	{
		if(this->stage==1)
			return false;
		else if(this->stage==0)
		{
			this->stage=1;
			double * newParameters = new double [this->NumParameters()];
			for(int i=0; i<3; i++)
			{
				newParameters[i]=this->Parameters[i];
			}
			for(int i=0; i< ::TunnelNumWells; i++)
			{
				newParameters[3*i+3]=this->Parameters[i+3];
				newParameters[3*i+4]=2.3;
				newParameters[3*i+5]=::DefaultWellLocations[i];
			}
			delete [] this->Parameters;
			this->Parameters=newParameters;
			return true;
		}
		else throw "error in TunnelwithGaussianInvParameters::Evolve : unexpected stage";
	}
	virtual bool Iterate(void)//for some non-continuous parameter, iterate through it. Intended to be used when it has never evolved
	{
		if(this->exp1==12 && this->exp2==10)
			return false;
		this->exp2+=2;
		if(this->exp2>=this->exp1)
		{
			this->exp2=4;
			this->exp1+=2;
		}
		return true;
	}
	virtual ParameterSet * clone() const
	{
		return new TunnelwithGaussianInvParameters(*this);
	}
	virtual void Write(std::ostream & out) const
	{
		out<<this->exp1<<" \t"<<this->exp2<<'\n';
		out<<this->stage<<'\n';
		out<<this->Parameters[0]<<'\n';
		out<<this->Parameters[1]<<" \t"<<this->Parameters[2]<<'\n';
		int end=(this->stage==1) ? 3+3*::TunnelNumWells : 3+::TunnelNumWells;
		for(int i=3; i<end; i++)
			out<<this->Parameters[i]<<" \t";
		out<<'\n';
	}
	virtual void Read(std::istream & in)
	{
		in>>this->exp1;
		in>>this->exp2;
		in>>this->stage;
		int end=(this->stage==1) ? 3+3*::TunnelNumWells : 3+::TunnelNumWells;
		delete [] this->Parameters;
		this->Parameters = new double [end];
		for(int i=0; i<end; i++)
			in>>this->Parameters[i];
	}
	virtual double RecommendedTolerence(void) const
	{
		return 1e-7;
	}
	virtual double RecommendedMaxEnergy(void) const//when optimizing specific volume, the smallest specific volume needs to be so small that this value is achieved
	{
		return 1e4;
	}
	virtual size_t NumParameters(void) const
	{
		return (this->stage==1) ? 3+3*::TunnelNumWells : 3+::TunnelNumWells;
	}
	virtual double MinDistance(void) const//when r<MinDistance, the potential should have no interesting behavior, guaranteed to fix unchanged during the life of this object
	{
		return 0.2;
	}
	virtual double MaxDistance(void) const//the Rcut of the potential should not exceed this, guaranteed to fix unchanged during the life of this object
	{
		return RCut;
	}
	virtual double SharpestRelativeSize(void)
	{
		return 0.02;
	}
};


#endif
