#include <cassert>

#include "MC_System.h"
#include "Potential.h"
#include "StructureOptimization.h"


#ifdef DEBUG
const double power1=0.04;
const double coeff1=20;
#else
const double power1=0.04;
const double coeff1=0.2;
#endif


bool ParticleSimulatedAnnealing::StopSignal=false;

////////////////////////////////////////////////////////////////////////////////////////
ParticleSimulatedAnnealing::ParticleSimulatedAnnealing(DimensionType Dimension, size_t NumParticle, double InitialDensity, int RandomSeed, double Pressure, Potential & pot, std::ostream * output, double CoolingScheduleRescale, double MinDistance, int AlgorithmLevel, bool UseSortedList) : Output(output), gen(RandomSeed), MinDistance(MinDistance)
{
	this->alpha=1.0- coeff1*std::pow( power1, Dimension)/NumParticle/CoolingScheduleRescale;
	if(alpha==1)
	{
		std::cerr<<"Error in ParticleSimulatedAnnealling::ParticleSimulatedAnnealing : System is too complex, unable to find a cooling schedule!";
		assert(false);
	}
	this->Tcut=0.00000100*std::pow(0.1, Dimension);

	this->Tinit=100000;
	this->CellMoveProbability=1/(5.0*NumParticle-4.0);

	this->Pressure=Pressure;
	size_t ExpectedStep=std::log(Tinit/Tcut)/std::log(alpha);

	//////////////////////
	//initialize this->pList
	if(NumParticle==0)
	{
		std::cerr<<"Error in MC_System::MC_System : NumParticle is 0\n";
		assert(false);
	}
	unsigned short cellrank = std::floor( std::pow(NumParticle/ParticlePerCell, static_cast<double>(1)/Dimension) );//about 1 particles per cell
	double initsize=std::pow(NumParticle/InitialDensity, static_cast<double>(1)/Dimension);
	std::vector<GeometryVector> base;
	for(DimensionType i=0; i<Dimension; i++)
	{
		base.push_back(GeometryVector(Dimension));
		base.back().x[i]=initsize;
		for(DimensionType j=i+1; j<Dimension; j++)
			base.back().x[j]=initsize*(gen.RandomDouble()-0.5);
	}
	this->pList=new Configuration(Dimension, &base[0], initsize/cellrank, UseSortedList);
	size_t trialcount=0;
	for(size_t i=0; i<NumParticle;)
	{
		GeometryVector temp;
		for(DimensionType j=0; j<Dimension; j++)
			temp.x[j]=gen.RandomDouble();
		bool TooNearNeighbor=false;
		this->pList->IterateThroughNeighbors(temp, MinDistance, [&TooNearNeighbor](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
		{
			TooNearNeighbor=true;
		}, &TooNearNeighbor);
		if(TooNearNeighbor==false)
		{
			this->pList->Insert("C", temp);
			i++;
		}
		trialcount++;
		if(trialcount%100000==99999)
		{
			GeometryVector tbase[ ::MaxDimension];
			for(DimensionType j=0; j<Dimension; j++)
				tbase[j]=this->pList->GetBasisVector(j)*2;
			this->pList->ChangeBasisVector(tbase);
		}
	}

	pot.SetConfiguration(*this->pList);
	double Enthalpy=pot.Energy()+Pressure*this->pList->PeriodicVolume();

	///////////////////////
	//initialize this->pSys
	if(AlgorithmLevel>1)
		this->pSys= new ParallelTempering<Configuration, Potential *> (*this->pList, gen.RandomDouble()*10000 , AlgorithmLevel, Enthalpy);
	else if(AlgorithmLevel==1)
		this->pSys= new MonteCarlo<Configuration, Potential *> (*this->pList, gen.RandomDouble()*10000, Enthalpy );
	else
		this->pSys=nullptr;

	if(this->Output!=nullptr)
	{
		(*this->Output).precision(14);
		(*this->Output)<<"Simulated Annealing Started, alpha="<<alpha<<", T_init="<<Tinit<<", T_cut="<<Tcut<<", Cell Move Probability="<<CellMoveProbability<<", Pressure="<<Pressure<<", Mindistance="<<this->MinDistance<<'\n';
	}
	this->RelaxStructureFunction= [] (Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance) ->void
	{
		RelaxStructure_ConjugateGradient(List, pot, Pressure, Switch, MinDistance, 10000);
	};
}
ParticleSimulatedAnnealing::ParticleSimulatedAnnealing(const Configuration & pInitConfig, int RandomSeed, double Pressure, Potential & pot, std::ostream * output, double CoolingScheduleRescale, double MinDistance, int AlgorithmLevel, bool UseSortedList) : Output(output), gen(RandomSeed), MinDistance(MinDistance)
{
	this->alpha=1.0- coeff1*std::pow( power1, pInitConfig.GetDimension())/pInitConfig.NumParticle()/CoolingScheduleRescale;
	if(alpha==1)
	{
		std::cerr<<"Error in ParticleSimulatedAnnealling::ParticleSimulatedAnnealing : System is too complex, unable to find a cooling schedule!";
		assert(false);
	}
	this->Tcut=0.00000100*std::pow(0.1, pInitConfig.GetDimension());

	this->Tinit=100000;
	this->CellMoveProbability=1/(5.0*pInitConfig.NumParticle()-4.0);

	this->Pressure=Pressure;
	size_t ExpectedStep=std::log(Tinit/Tcut)/std::log(alpha);

	//////////////////////
	//initialize this->pList
	this->pList=new Configuration(pInitConfig);
	pot.SetConfiguration(*this->pList);
	double Enthalpy=pot.Energy()+Pressure*this->pList->PeriodicVolume();

	///////////////////////
	//initialize this->pSys
	if(AlgorithmLevel>1)
		this->pSys= new ParallelTempering<Configuration, Potential *> (*this->pList, gen.RandomDouble()*10000 , AlgorithmLevel, Enthalpy);
	else if(AlgorithmLevel==1)
		this->pSys= new MonteCarlo<Configuration, Potential *> (*this->pList, gen.RandomDouble()*10000, Enthalpy );
	else
		this->pSys=nullptr;

	this->RelaxStructureFunction= [] (Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance) ->void
	{
		RelaxStructure_ConjugateGradient(List, pot, Pressure, Switch, MinDistance, 10000);
	};
}

class StopMCAnnealing
{
public:
	StopMCAnnealing()
	{}
};
//const double ZeroMoveStopTolerence = 1e-10;
void ParticleSimulatedAnnealing::Anneal(Potential & pot)
{
	if(this->Output!=nullptr)
	{
		(*this->Output).precision(14);
		(*this->Output)<<"Simulated Annealing Started, alpha="<<alpha<<", T_init="<<Tinit<<", T_cut="<<Tcut<<", Cell Move Probability="<<CellMoveProbability<<", Pressure="<<Pressure<<", Mindistance="<<this->MinDistance<<'\n';
	}
	AtomMove ma;
	CellMove mb(this->Pressure);
	ma.MinDistance=this->MinDistance;
	mb.MinDistance=this->MinDistance;

	//debug temp
	//this->CellMoveProbability=1.0;
	CombinedMove<Configuration, Potential *> comb(ma, mb, 1.0-this->CellMoveProbability, this->gen.RandomDouble()*10000);

	double p = this->Pressure;
	std::ostream * op= Output;
	double prevT;
	bool RelaxCell=(this->CellMoveProbability!=0.0);
	if(this->pSys != nullptr)
	{
		pot.SetConfiguration(*this->pList);
		this->pSys->RecordEnergy=pot.Energy()+p*this->pList->PeriodicVolume();

		ThermoCool cool2(this->alpha, this->Tinit, this->Tcut, 0);
		//ExpCool cool2(this->alpha, this->Tinit, this->Tcut, 0);
		size_t * pMoveCount = & this->pSys->MoveCount;
		double prevCheckRelaxStructureT=this->Tinit;
		double prevCheckRelaxStructureH=this->pSys->RecordEnergy;
		double MinDist=this->MinDistance;
		size_t GotSameLocalCount=0;
		bool * pStopSignal=& ParticleSimulatedAnnealing::StopSignal;
		auto * pRelaxFunction = & this->RelaxStructureFunction;

		try
		{
			this->pSys->Anneal(&pot, cool2, comb, [&](double & Energy, double Temperature, Configuration & l)
			{
				pot.SetConfiguration(l);
				double newE=pot.Energy()+p*l.PeriodicVolume();
				if(newE<1e5 && std::abs((newE-Energy)/(Energy)) > 1e-8 && std::abs((newE-Energy)) > 1e-8)
					std::cerr<<"Error in MC_System.cpp: Energy not consistent, newE="<<newE<<", record E="<<Energy<<'\n';
				Energy=newE;
				if(Temperature!=prevT)
				{
					if(op!=nullptr)
						(*op)<<*pMoveCount<<" \t"<<Temperature<<" \t"<<l.PeriodicVolume()<<" \t"<<Energy<<'\n';
					prevT=Temperature;
				}
				if((*(pStopSignal))==true)
					throw StopMCAnnealing();

				if(Temperature<0.5*prevCheckRelaxStructureT || (Temperature<0.85*prevCheckRelaxStructureT&&GotSameLocalCount>0))
				{
					//debug temp
					//std::string prefix("Temperature");
					//prefix+=std::to_string( (long double)(Temperature) );
					//::Plot(prefix, l);
					//check if we get the same enthalpy after Relaxing
					Configuration temp(l);
					if(RelaxCell)
						(*pRelaxFunction)(temp, pot, p, 2, MinDist);
					else
						(*pRelaxFunction)(temp, pot, p, 0, MinDist);
					pot.SetConfiguration(temp);
					double tempEnthalpy=pot.Energy()+p*temp.PeriodicVolume();
					if(op!=nullptr)
						(*op)<<"Prev local minimum H="<<prevCheckRelaxStructureH<<", Now local minimum H="<<tempEnthalpy;
					if( std::abs((tempEnthalpy-prevCheckRelaxStructureH)/prevCheckRelaxStructureH)<1e-7 || std::abs((tempEnthalpy-prevCheckRelaxStructureH))<1e-7 )//we get the same local minima as previous
					{
						GotSameLocalCount++;
						if(op!=nullptr)
							(*op)<<", Achieved the same value"<<GotSameLocalCount<<"times, ";
						if(GotSameLocalCount > 6)
						{
							if(op!=nullptr)
								(*op)<<"Stop Cooling\n";
							throw StopMCAnnealing();
						}
						else
							if(op!=nullptr)
								(*op)<<"Continue Cooling\n";

						prevCheckRelaxStructureT=Temperature;
					}
					else
					{
						if(op!=nullptr)
							(*op)<<", Continue Cooling\n";
						GotSameLocalCount=0;
						prevCheckRelaxStructureH=tempEnthalpy;
						prevCheckRelaxStructureT=Temperature;
					}
					if(op!=nullptr)
						(*op).flush();
				}
				pot.SetConfiguration(l);
			});
		}
		catch(StopMCAnnealing & a)
		{}
	}

	if(ParticleSimulatedAnnealing::StopSignal==false)
	{
		pot.SetConfiguration(*this->pList);
		if(op!=nullptr)
		{
			(*op)<<"Optimization Start, Starting size:"<<this->pList->PeriodicVolume()<<", Enthalpy:"<<pot.Energy()+this->Pressure*this->pList->PeriodicVolume()<<'\n';
		}
		ParallelTempering<Configuration, Potential *> * pPT = dynamic_cast<ParallelTempering<Configuration, Potential *> *>(this->pSys);
		if(pPT==nullptr)
		{
			if(RelaxCell)
				this->RelaxStructureFunction(*this->pList, pot, p, 2, this->MinDistance);
			else
				this->RelaxStructureFunction(*this->pList, pot, p, 0, this->MinDistance);
		}
		else
		{
			pot.SetConfiguration(*this->pList);
			double Enthalpy=pot.Energy()+this->Pressure*this->pList->PeriodicVolume();
			for(auto iter=pPT->vMonteCarlos.begin(); iter!=pPT->vMonteCarlos.end(); iter++)
			{
				Configuration now(*iter->pSys);
				if(RelaxCell)
					this->RelaxStructureFunction(now, pot, p, 2, this->MinDistance);
				else
					this->RelaxStructureFunction(now, pot, p, 0, this->MinDistance);
				pot.SetConfiguration(now);
				double nowEnthalpy=pot.Energy()+this->Pressure*now.PeriodicVolume();
				if(nowEnthalpy<Enthalpy)
				{
					Enthalpy=nowEnthalpy;
					(*this->pList)=now;
				}
			}
		}
		pot.SetConfiguration(*this->pList);
		if(op!=nullptr)
		{
			(*op)<<"Optimization Finish, Ending size:"<<this->pList->PeriodicVolume()<<", Enthalpy:"<<pot.Energy()+this->Pressure*this->pList->PeriodicVolume()<<'\n';
		}
	}
}

