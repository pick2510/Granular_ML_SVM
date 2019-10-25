#ifndef MCMOVE_INCLUDED
#define MCMOVE_INCLUDED

#include "etc.h"
#include "PeriodicCellList.h"
#include "Potential.h"
#include "CollectiveCoordinatePotential.h"
#include "RandomGenerator.h"
#include <iostream>
#include <cassert>
#include <exception>
const size_t AcceptanceSampleSize=1000;
const double MaxAcceptance=0.7;
const double MinAcceptance=0.3;


const double MaxAspectRatioOfCell=1.9;
const double MaxCellSideLengthSquared=10000.0;


template <typename System, typename AdditionalInfo>
class MCMove
{
public:
	//Generate a random move, return Energy difference. the PreFactor is also set if it is not 1.
	virtual double DeltaEnergy(System & sys, AdditionalInfo add, RandomGenerator & gen, double & PreFactor, bool LockStepSize, double CurrentEnergy) =0;
	virtual void Accept(System & sys, AdditionalInfo add) =0;
	virtual MCMove<System, AdditionalInfo> * clone() =0;
};

class AtomMove : public MCMove<Configuration, Potential *>
{
private:
	size_t TrialCount, AcceptCount;

	size_t moved;
	GeometryVector newcoord;

public:
	//the most recently calculated acceptance ratio, updated by the class
	double AcceptRatio;
	double sigma;
	double MyMinAcceptance;

	double MinDistance;//abandon this move immediately if particles are within this distance
	AtomMove()
	{
		this->sigma=0.1;
		this->TrialCount=0;
		this->AcceptCount=0;
		this->MyMinAcceptance= ::MinAcceptance;
		this->MinDistance=0;
	}

	//calculate the delta E for a given move
	//move is : move the particle (this->moved) so that its new relative coordinate is newcoord
	double CalculateDeltaEnergy(Configuration & list, Potential * pot, const GeometryVector & newcoord, double CurrentEnergy)
	{
		PairPotential * pPair = dynamic_cast<PairPotential *>(pot);
		forMCPotential * pForMCPotential = dynamic_cast<forMCPotential *>(pot);
		CombinedPotential * pComb = dynamic_cast<CombinedPotential *>(pot);
		//if(pCC!=nullptr && pForMCPotential==nullptr)
		//	std::cout<<"Warning in AtomMove : ShiftedCCPotential does not have good performance for this purpose. Use forMCPotential instead!\n";
		if(pPair!=nullptr)
		{
			double deltaEnergy=0;
			AtomInfo a1 = list.GetCharacteristics(moved);
			list.IterateThroughNeighbors(newcoord, pPair->Rcut, [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) -> void 
			{
				AtomInfo a2 = list.GetCharacteristics(SourceAtom);
				if(SourceAtom!=moved)
					deltaEnergy+=pPair->EnergyFormula(shift, a1, a2);
			});
			list.IterateThroughNeighbors(moved, pPair->Rcut, [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) -> void 
			{
				AtomInfo a2 = list.GetCharacteristics(SourceAtom);
				if(SourceAtom!=moved)
					deltaEnergy -= pPair->EnergyFormula(shift, a1, a2);
			});
			return deltaEnergy;
		}
		else if (pForMCPotential!=nullptr)
		{
			GeometryVector prevCart=list.GetCartesianCoordinates(this->moved);
			GeometryVector afterCart=list.RelativeCoord2CartesianCoord(this->newcoord);
			return pForMCPotential->TryMove(prevCart, afterCart);
		}
		else if (pComb!=nullptr)
		{
			return this->CalculateDeltaEnergy(list, pComb->p1, newcoord, CurrentEnergy)+this->CalculateDeltaEnergy(list, pComb->p2, newcoord, CurrentEnergy);
		}
		else
		{
			GeometryVector backup=list.GetRelativeCoordinates(this->moved);
			pot->SetConfiguration(list);
			double E1=pot->Energy();

			list.MoveParticle(this->moved, this->newcoord); 
			pot->SetConfiguration(list);
			double E2=pot->Energy();
			list.MoveParticle(this->moved, backup); 

			return E2-E1;
		}
	}
	virtual double DeltaEnergy(Configuration & list, Potential * pot, RandomGenerator & gen, double & PreFactor, bool LockStepSize, double CurrentEnergy)//Generate a random move, return Energy difference
	{
		if(this->TrialCount>=::AcceptanceSampleSize )
		{
			//TODO : THIS MAY BREAK DETAILED BALANCE. FIND A WAY TO FIX IT.
			AcceptRatio=static_cast<double>(this->AcceptCount)/this->TrialCount;//acceptance rate
			if(LockStepSize==false)
			{
				if(AcceptRatio> ::MaxAcceptance && this->sigma<0.5)
				{
					this->sigma*=1.1;
					//std::cout<<"In AtomMove, step size changed to"<<this->sigma<<"\n";
				}
				else if(AcceptRatio < this->MyMinAcceptance)
				{
					this->sigma*=0.9;
					//std::cout<<"In AtomMove, step size changed to"<<this->sigma<<"\n";
				}
				//else
					//std::cout<<"In AtomMove, keep step size"<<this->sigma<<"\n";
			}
			//debug temp
			//std::cout<<"Sigma is:"<<this->sigma<<", Accept Count is:"<<this->AcceptCount<<'\n';
			this->TrialCount=0;
			this->AcceptCount=0;
		}
		this->TrialCount++;

		bool TooNearNeighborFound=true;

		size_t ReTrialCount=0;
		this->moved=list.GetRandomParticle(gen);
		this->newcoord=list.GetRelativeCoordinates(this->moved);
		DimensionType dim=list.GetDimension();
		for(int i=0; i<dim; i++)
			this->newcoord.x[i]+=this->sigma*(gen.RandomDouble()-0.5);

		TooNearNeighborFound=false;
		auto moved=this->moved;

		if (MinDistance > 0.0)
		{
			list.IterateThroughNeighbors(newcoord, this->MinDistance, [&TooNearNeighborFound, &moved](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) -> void
			{
				if (SourceAtom != moved)
					TooNearNeighborFound = true;
			}, &TooNearNeighborFound);
		}

		//new_idea
		if(TooNearNeighborFound)
			return ::MaxEnergy;
		//restore?ReTrialCount++;

		return CalculateDeltaEnergy(list, pot, newcoord, CurrentEnergy);
	}

	//check if the potential needs specific operation for accepting the move
	void AcceptCheck(Potential * pot)
	{
		forMCPotential * pForMCPotential = dynamic_cast<forMCPotential *>(pot);
		CombinedPotential * pComb= dynamic_cast<CombinedPotential *>(pot);
		if(pForMCPotential!=nullptr)
			pForMCPotential->AcceptMove();
		else if(pComb!=nullptr)
		{
			this->AcceptCheck(pComb->p1);
			this->AcceptCheck(pComb->p2);
		}
	}
	virtual void Accept(Configuration & list, Potential * pot)
	{
		this->AcceptCount++;

		list.MoveParticle(this->moved, this->newcoord); 

		this->AcceptCheck(pot);
	}
	virtual MCMove<Configuration, Potential *> * clone()
	{
		return new AtomMove(*this);
	}
};

class CellMove : public MCMove<Configuration, Potential *>
{
private:
	size_t TrialCount, AcceptCount;
	double Pressure;

	GeometryVector newvec[::MaxDimension];
public:
	double sigma;
	double MinDistance;//abandon this move immediately if particles are within this distance
	CellMove()
	{
		this->sigma=3e-2;
		this->TrialCount=0;
		this->AcceptCount=0;
		this->Pressure=0;
		this->MinDistance=0;
	}
	CellMove(double pressure)
	{
		this->sigma=3e-2;
		this->TrialCount=0;
		this->AcceptCount=0;
		this->Pressure=pressure;
		this->MinDistance=0;
	}
	virtual double DeltaEnergy(Configuration & list, Potential * pot, RandomGenerator & gen, double & PreFactor, bool LockStepSize, double CurrentEnergy)//Generate a random move, return Energy difference
	{
		if(this->TrialCount>=::AcceptanceSampleSize && LockStepSize==false)
		{
			double arate=static_cast<double>(this->AcceptCount)/this->TrialCount;//acceptance rate
			//std::cout<<"arate is:"<<arate<<'\n';
			if(arate> ::MaxAcceptance)
			{
				//std::cout<<"Warning in CellMove : step size changed!\n";
				this->sigma*=1.1;
				if(this->sigma>0.3)
					this->sigma=0.3;
				//std::cout<<"Sigma is:"<<this->sigma<<'\n';
			}
			else if(arate < ::MinAcceptance)
			{
				//std::cout<<"Warning in CellMove : step size changed!\n";
				this->sigma*=0.9;
				//std::cout<<"Sigma is:"<<this->sigma<<'\n';
			}
			this->TrialCount=0;
			this->AcceptCount=0;
			//std::cout<<"Sigma is:"<<this->sigma<<'\n';
		}
		if(this->TrialCount==::AcceptanceSampleSize/10)
		{
			double arate=static_cast<double>(this->AcceptCount)/this->TrialCount;//acceptance rate
			//std::cout<<"arate is:"<<arate<<'\n';
			if(arate < ::MinAcceptance/2)
			{
				this->sigma*=0.5;
				this->TrialCount=0;
				this->AcceptCount=0;
				//std::cout<<"Sigma is:"<<this->sigma<<'\n';
			}
		}
		this->TrialCount++;

		auto dim=list.GetDimension();
		bool CellWeired;
		bool TooNearNeighborFound=true;
		Configuration newlist(list);
		size_t ReTrialTime =0;//don't try the following loop too much, it's not terrible even if we don't have a perfect move
		//pot->SetConfiguration(list);
		//double OldH = pot->Energy() + list.PeriodicVolume()*this->Pressure;
		double OldH = CurrentEnergy;
		do
		{
			CellWeired=false;
			double L2Max=0, L2Min=1000000000;
			for(DimensionType i=0; i<dim; i++)
				newvec[i]=list.GetBasisVector(i);
			for(DimensionType i=0; i<dim; i++)
			{
				for(int j=0; j<list.GetDimension(); j++)
				{
					double length=std::sqrt(newvec[i].Modulus2());
                    newvec[i].x[j]+=this->sigma*(gen.RandomDouble()-0.5)*(length);
				}
				double dtemp=newvec[i].Modulus2();
				if(dtemp>L2Max)
					L2Max=dtemp;
				if(dtemp<L2Min)
					L2Min=dtemp;
			}
			double OrigVol=list.PeriodicVolume();
			double NewVol=::Volume(newvec, dim);
			if(std::abs((NewVol-OrigVol)/OrigVol) > this->sigma )
				CellWeired=true;
			if( (L2Max/L2Min)> (MaxAspectRatioOfCell)*(MaxAspectRatioOfCell) )
				CellWeired=true;
			if(L2Max > MaxCellSideLengthSquared)
				CellWeired=true;

			//search for too near neighbors
			newlist.ChangeBasisVector(newvec);
			TooNearNeighborFound=false;
			if (MinDistance > 0.0)
			{
				for (size_t i = 0; i < newlist.NumParticle(); i++)
				{
					//Configuration::particle * pa=newlist.GetParticle(i);
					newlist.IterateThroughNeighbors(i, this->MinDistance, [&TooNearNeighborFound, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) -> void
					{
						if (SourceAtom != i || LatticeShift.Modulus2() != 0)
							TooNearNeighborFound = true;
					}, &TooNearNeighborFound);
				}
			}
			//new_idea
			ReTrialTime++;
			if( (TooNearNeighborFound||CellWeired) && ReTrialTime>100 )
				return ::MaxEnergy;
		}
		while( (CellWeired || TooNearNeighborFound) );

		double NewH;
		try
		{
			pot->SetConfiguration(newlist);
			NewH = pot->Energy() + newlist.PeriodicVolume()*this->Pressure;
		}
		catch(std::bad_alloc & a)
		{
			std::cerr<<"Warning : not enough memory in cell resize, abandon move!";
			return ::MaxEnergy;
		}

		PreFactor=std::pow(newlist.PeriodicVolume()/list.PeriodicVolume(), (double)(list.NumParticle()));

		double DeltaE= NewH-OldH;
		//std::cout <<" \t"<< PreFactor << " \t" << DeltaE << '\n';
		return DeltaE;
	}
	virtual void Accept(Configuration & list, Potential * pot)
	{
		list.ChangeBasisVector(newvec);
		list.TryRefineBasisVectors();
		this->AcceptCount++;
	}
	virtual MCMove<Configuration, Potential *> * clone()
	{
		return new CellMove(*this);
	}
};


template <typename System, typename AdditionalInfo>
class CombinedMove : public MCMove<System, AdditionalInfo>
{
	MCMove<System, AdditionalInfo> * m1;
	MCMove<System, AdditionalInfo> * m2;
	double Move1Probability;
	RandomGenerator gen;
	bool IsMoving1;

public:
	CombinedMove(MCMove<System, AdditionalInfo> & m1, MCMove<System, AdditionalInfo> & m2, double Move1Probability, int RandomSeed) : m1(m1.clone()), m2(m2.clone()), Move1Probability(Move1Probability), gen(RandomSeed)
	{
		assert(Move1Probability>=0);
		assert(Move1Probability<=1);
	}
	~CombinedMove()
	{
		delete this->m1;
		delete this->m2;
	}
	virtual double DeltaEnergy(System & sys, AdditionalInfo add, RandomGenerator & gen, double & PreFactor, bool LockStepSize, double CurrentEnergy)//Generate a random move, return Energy difference
	{
		if(gen.RandomDouble()<this->Move1Probability)
		{
			this->IsMoving1=true;
			return m1->DeltaEnergy(sys, add, gen, PreFactor, LockStepSize, CurrentEnergy);
		}
		else
		{
			this->IsMoving1=false;
			return m2->DeltaEnergy(sys, add, gen, PreFactor, LockStepSize, CurrentEnergy);
		}
	}
	virtual void Accept(System & sys, AdditionalInfo add)
	{
		if(this->IsMoving1==true)
			this->m1->Accept(sys, add);
		else
			this->m2->Accept(sys, add);
	}
	virtual MCMove<System, AdditionalInfo> * clone()
	{
		return new CombinedMove(*this->m1, *this->m2, this->Move1Probability, this->gen.RandomDouble()*10000);
	}
};


#endif
