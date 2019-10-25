#include "MD_System.h"
#include "Plots.h"
#include "StructureOptimization.h"
#include <vector>
#include <cmath>
#include "etc.h"


//debug temp
#include <float.h>


void ParticleMolecularDynamics::InitNoseHoover(double Q)
{
	this->MassQ=Q;
	this->xi=0;
}

ParticleMolecularDynamics::ParticleMolecularDynamics(const Configuration & Position, double TimeStep, double Mass)
	: Position(Position)
	, Velocity(Position.NumParticle(), GeometryVector(Position.GetDimension()))
	, Acceleration(Position.NumParticle(), GeometryVector(Position.GetDimension()))
	, Mass(Position.NumParticle(), Mass)
{
	this->NumParticle=Position.NumParticle();
	this->Dimension=Position.GetDimension();
	this->TimeStep=TimeStep;
	this->Time=0;
	this->MassQ=0;
}

ParticleMolecularDynamics::~ParticleMolecularDynamics()
{
}

ParticleMolecularDynamics::ParticleMolecularDynamics(std::istream & ifile) : Position(ifile)
{
	ifile.read( (char*)(&TimeStep), sizeof(TimeStep) );
	ifile.read( (char*)(&Time), sizeof(Time) );
	ifile.read( (char*)(&NumParticle), sizeof(NumParticle) );
	ifile.read( (char*)(&Dimension), sizeof(Dimension) );
	ifile.read( (char*)(&xi), sizeof(xi) );
	ifile.read( (char*)(&MassQ), sizeof(MassQ) );
	this->Mass.resize(this->NumParticle);
	ifile.read( (char*)(&Mass[0]), sizeof(Mass[0])*this->NumParticle );
	this->Velocity.resize(this->NumParticle);
	ifile.read( (char*)(&Velocity[0]), sizeof(Velocity[0])*this->NumParticle );
	this->Acceleration.resize(this->NumParticle);
	ifile.read( (char*)(&Acceleration[0]), sizeof(Acceleration[0])*this->NumParticle );
}
void ParticleMolecularDynamics::WriteBinary(std::ostream & ofile)
{
	this->Position.WriteBinary(ofile);
	ofile.write( (char*)(&TimeStep), sizeof(TimeStep) );
	ofile.write( (char*)(&Time), sizeof(Time) );
	ofile.write( (char*)(&NumParticle), sizeof(NumParticle) );
	ofile.write( (char*)(&Dimension), sizeof(Dimension) );
	ofile.write( (char*)(&xi), sizeof(xi) );
	ofile.write( (char*)(&MassQ), sizeof(MassQ) );
	ofile.write( (char*)(&Mass[0]), sizeof(Mass[0])*this->NumParticle );
	ofile.write( (char*)(&Velocity[0]), sizeof(Velocity[0])*this->NumParticle );
	ofile.write( (char*)(&Acceleration[0]), sizeof(Acceleration[0])*this->NumParticle );
}


void ParticleMolecularDynamics::_GetAcc(Potential & potential)
{
	//potential.Force(this->Position, this->Velocity, this->Acceleration->GetParticle(i), i);
	//for(int j=0; j<this->Dimension; j++)
	//	this->Acceleration->GetParticle(i)[j]/=this->Mass;
	potential.AllForce(this->Acceleration);
	for(int i=0; i<this->Acceleration.size(); i++)
		this->Acceleration[i].MultiplyFrom( 1.0/this->Mass[i] );
}

double ParticleMolecularDynamics::GetPotentialEnergy(Potential & potential)
{
	potential.SetConfiguration(this->Position);
	return potential.Energy();
}	

//TODO : Why does this function break equilibrium and greatly increase the potential energy?
void ParticleMolecularDynamics::FindOptimalTimeStep(Potential & potential, double & TimeStep, double TargetDeltaForce, double kT)
{
	RandomGenerator gen(1111);
	for(;;)
	{
		// find MaxDeltaForce
		this->MaxDeltaForce=0.0;
		this->TimeStep=TimeStep;
		for(int i=0; i<100; i++)
		{
			this->SetRandomSpeed(kT, gen);
			this->Evolve(10, potential);
			//std::cout<<"E_k_Ev="<<this->GetKineticEnergy();
			this->Evolve_RecordDeltaForce(10, potential);
			//std::cout<<"E_k_EvR="<<this->GetKineticEnergy()<<'\n';
		}
		if( this->MaxDeltaForce>5*TargetDeltaForce )
			TimeStep/=2;
		else if(this->MaxDeltaForce>2*TargetDeltaForce)
			TimeStep/=1.2;
		else if(this->MaxDeltaForce>1.2*TargetDeltaForce)
			TimeStep/=1.05;
		else if(this->MaxDeltaForce<0.2*TargetDeltaForce)
			TimeStep*=2;
		else if(this->MaxDeltaForce<0.5*TargetDeltaForce)
			TimeStep*=1.2;
		else if(this->MaxDeltaForce<0.8*TargetDeltaForce)
			TimeStep*=1.05;
		else
			break;
		std::cout<<"Max Delta force="<<this->MaxDeltaForce<<", TimeStep adjusted to:"<<TimeStep<<'\n';
		logfile<<"Max Delta force="<<this->MaxDeltaForce<<", TimeStep adjusted to:"<<TimeStep<<'\n';
		if(this->MaxDeltaForce>1e10)
		{
			std::cout<<"max delta force too large, relax structure!\n";
			logfile<<"max delta force too large, relax structure!\n";
			RelaxStructure_NLOPT(this->Position, potential, 0.0, 0, 0.0);
		}
	}
	this->TimeStep=TimeStep;
}

void ParticleMolecularDynamics::Evolve(size_t repeat, Potential & potential)//flag=0: serial, flag=1:omp parallel;
{
//#pragma omp parallel default(shared)
	{
		for(size_t u=0; u<repeat; u++)
		{
			//calculate v(t+0.5dt)
//#pragma omp for
			for(long i=0; i<this->NumParticle; i++)
				this->Velocity[i].AddFrom(0.5*this->TimeStep*this->Acceleration[i]);

			//calculate x(t+dt)
//#pragma omp single
			for(long i=0; i<this->NumParticle; i++)
			{
				//this->Position->GetParticle(i)[j]+=this->TimeStep*this->Velocity->GetParticle(i)[j];
				GeometryVector newCart=this->Position.GetCartesianCoordinates(i)+this->TimeStep*this->Velocity[i];
				this->Position.MoveParticle( i, this->Position.CartesianCoord2RelativeCoord(newCart) );
			}
			//calculate a(t+dt)

//#pragma omp single
			{
				potential.SetConfiguration(this->Position);
				this->_GetAcc(potential);
			}

			//calculate v(t+dt)
//#pragma omp for
			for(long i=0; i<this->NumParticle; i++)
				this->Velocity[i].AddFrom( 0.5*this->TimeStep*this->Acceleration[i] );

			//add time
//#pragma omp single
			this->Time+=this->TimeStep;
		}
	}
}

void ParticleMolecularDynamics::Evolve_AutoTimeStep(size_t repeat, Potential & potential, double DeltaLogEnergyPerStep)
{
	double Ebefore=this->GetKineticEnergy()+this->GetPotentialEnergy(potential);
	this->Evolve(repeat, potential);
	double Eafter=this->GetKineticEnergy()+this->GetPotentialEnergy(potential);
	double dlogEperStep= std::abs( std::log(Ebefore/Eafter)/repeat );

	if (dlogEperStep>5 * DeltaLogEnergyPerStep)
	{
		this->TimeStep *= 0.5;
	}
	else if (dlogEperStep>3 * DeltaLogEnergyPerStep)
	{
		this->TimeStep *= 0.9;
	}
	else if (dlogEperStep>2 * DeltaLogEnergyPerStep)
	{
		this->TimeStep*=0.95;
	}
	else if(dlogEperStep<0.01*DeltaLogEnergyPerStep)
	{
		this->TimeStep*=1.20;
	}
	else if(dlogEperStep<0.5*DeltaLogEnergyPerStep)
	{
		this->TimeStep*=1.05;
	}

	//debug temp
	//std::cout<<"TimeStep Adjusted to "<<this->TimeStep<<'\n';

}

void ParticleMolecularDynamics::AndersonEvolve(size_t repeat, Potential & potential, double kT, double ResetSpeedProbability, RandomGenerator & gen)
{
	double ResetCount=0.0;
//#pragma omp parallel default(shared)
	{
		for(size_t u=0; u<repeat; u++)
		{
			//calculate v(t+0.5dt)
//#pragma omp for
			for(long i=0; i<this->NumParticle; i++)
				this->Velocity[i].AddFrom(0.5*this->TimeStep*this->Acceleration[i]);

			//calculate x(t+dt)
//#pragma omp single
			for(long i=0; i<this->NumParticle; i++)
			{
				//this->Position->GetParticle(i)[j]+=this->TimeStep*this->Velocity->GetParticle(i)[j];
				GeometryVector newCart=this->Position.GetCartesianCoordinates(i)+this->TimeStep*this->Velocity[i];
				this->Position.MoveParticle( i, this->Position.CartesianCoord2RelativeCoord(newCart) );
			}
			//calculate a(t+dt)
//#pragma omp single
			{
				potential.SetConfiguration(this->Position);
				this->_GetAcc(potential);
			}

			//calculate v(t+dt)
//#pragma omp for
			for(long i=0; i<this->NumParticle; i++)
				this->Velocity[i].AddFrom( 0.5*this->TimeStep*this->Acceleration[i] );

			//add time
//#pragma omp single
			{
				this->Time+=this->TimeStep;
				ResetCount+=ResetSpeedProbability;
				for(; ResetCount>0; ResetCount-=1.0)
				{
					size_t i=std::floor(gen.RandomDouble()*this->NumParticle);
					double sigma=std::sqrt(kT/this->Mass[i]);
					for (DimensionType j = 0; j < this->Dimension; j++)
					{
						double temp;
						do
						{
							temp = gen.RandomDouble();
						} while (temp <= 0);
						this->Velocity[i].x[j] = std::sqrt(-2 * std::log(temp))*std::sin(2 * pi*gen.RandomDouble())*sigma;
					}
				}
				//if(gen.RandomDouble() > ResetSpeedProbability)
				//{
				//	size_t i=std::floor(gen.RandomDouble()*this->NumParticle);
				//	double sigma=std::sqrt(kT/this->Mass[i]);
				//	for(DimensionType j=0; j<this->Dimension; j++)
				//		this->Velocity[i].x[j]=std::sqrt(-2*std::log(gen.RandomDouble()))*std::sin(2*pi*gen.RandomDouble())*sigma;
				//}
			}
		}
	}
}


void ParticleMolecularDynamics::Evolve_RecordDeltaForce(size_t repeat, Potential & potential)//flag=0: serial, flag=1:omp parallel;
{
	{
		for(size_t u=0; u<repeat; u++)
		{
			//calculate v(t+0.5dt)
			for(long i=0; i<this->NumParticle; i++)
				this->Velocity[i].AddFrom(0.5*this->TimeStep*this->Acceleration[i]);

			//calculate x(t+dt)
			for(long i=0; i<this->NumParticle; i++)
			{
				//this->Position->GetParticle(i)[j]+=this->TimeStep*this->Velocity->GetParticle(i)[j];
				GeometryVector newCart=this->Position.GetCartesianCoordinates(i)+this->TimeStep*this->Velocity[i];
				this->Position.MoveParticle( i, this->Position.CartesianCoord2RelativeCoord(newCart) );
			}
			//calculate a(t+dt)
			potential.SetConfiguration(this->Position);
			//double f0a=this->Acceleration[0].x[0];
			std::vector<GeometryVector> fa = this->Acceleration;
			this->_GetAcc(potential);

			double AverageForce=0.0;
			for(size_t i=0; i<this->NumParticle; i++)
				for(DimensionType j=0; j<this->Dimension; j++)
					AverageForce+=std::abs(fa[i].x[j]);
			AverageForce=AverageForce/this->NumParticle/this->Dimension;
			for(size_t i=0; i<this->NumParticle; i++)
			{
				for(DimensionType j=0; j<this->Dimension; j++)
				{
					double f0a=fa[i].x[j];
					double f0b=this->Acceleration[i].x[j];
					double DeltaForcePerStep=std::abs((f0a-f0b)/AverageForce);
					if(this->MaxDeltaForce<DeltaForcePerStep)
						this->MaxDeltaForce=DeltaForcePerStep;
				}
			}

			//calculate v(t+dt)
			for(long i=0; i<this->NumParticle; i++)
				this->Velocity[i].AddFrom( 0.5*this->TimeStep*this->Acceleration[i] );

			//add time
			this->Time+=this->TimeStep;
		}
	}
}
void ParticleMolecularDynamics::NoseHooverEvolve(size_t repeat, Potential & potential, double kT)//flag=0: serial, flag=1:omp parallel;
{
	//SnapShot * halfvelocity, * newvelocity;
	double oldxi, temp;
	std::vector<GeometryVector> halfvelocity, newvelocity(this->NumParticle, GeometryVector(this->Dimension));
//#pragma omp parallel default(shared)
	{
		for(size_t u=0; u<repeat; u++)
		{
//#pragma omp master
			halfvelocity = std::vector<GeometryVector>(this->NumParticle, GeometryVector(this->Dimension));
//#pragma omp barrier
			//calculate v(t+0.5dt)
//#pragma omp for
			for(long i=0; i<this->NumParticle; i++)
				this->Velocity[i].AddFrom(0.5*this->TimeStep*this->Acceleration[i]);

			//calculate x(t+dt)
//#pragma omp single
			for(long i=0; i<this->NumParticle; i++)
			{
				//this->Position->GetParticle(i)[j]+=this->TimeStep*this->Velocity->GetParticle(i)[j];
				GeometryVector newCart=this->Position.GetCartesianCoordinates(i)+this->TimeStep*this->Velocity[i];
				this->Position.MoveParticle( i, this->Position.CartesianCoord2RelativeCoord(newCart) );
			}
			//calculate a(t+dt)

//#pragma omp single
			{
				potential.SetConfiguration(this->Position);
				this->_GetAcc(potential);
			}

			//calculate v(t+dt)
//#pragma omp single
			{
				//newvelocity = new SnapShot(this->Dimension, this->NumParticle, this->Velocity->SideLength);
				oldxi=this->xi;
				for(int m=0; m<20; m++)
				{
					for(long i=0; i<this->NumParticle; i++)
						newvelocity[i]=(halfvelocity[i]+this->Acceleration[i]*0.5*this->TimeStep)*(1.0/(1+this->xi*0.5*this->TimeStep));

					temp=0+oldxi;
					for(long i=0; i<this->NumParticle; i++)
						temp+=(this->Mass[i]*this->Velocity[i].Modulus2()+this->Mass[i]*newvelocity[i].Modulus2()-2*kT)/this->MassQ*0.5*this->TimeStep;
					this->xi=temp;
				}
				//this->Velocity=newvelocity;
				std::swap(this->Velocity, newvelocity);

				//add time
				this->Time+=this->TimeStep;
			}
		}
	}
}
double ParticleMolecularDynamics::GetKineticEnergy(void)
{
	double result=0;
	for(size_t i=0; i<this->NumParticle; i++)
		result+=this->Velocity[i].Modulus2()*this->Mass[i];
	result*=0.5;
	return result;
}

void ParticleMolecularDynamics::SetRandomSpeed(double kT, RandomGenerator & gen)
{
	for(size_t i=0; i<this->NumParticle; i++)
	{
		double sigma=std::sqrt(kT/this->Mass[i]);
		for (DimensionType j = 0; j < this->Dimension; j++)
		{
			this->Velocity[i].x[j] = std::sqrt(-2 * std::log(gen.RandomDouble()))*std::sin(2 * pi*gen.RandomDouble())*sigma;

			//gen.RandomDouble could return exactly zero, so check for blowing up
			while(std::abs(this->Velocity[i].x[j])>1e30)
				this->Velocity[i].x[j] = std::sqrt(-2 * std::log(gen.RandomDouble()))*std::sin(2 * pi*gen.RandomDouble())*sigma;
		}
	}
	return;
}
void ParticleMolecularDynamics::RemoveOverallMomentum()
{
	GeometryVector op;

	for(size_t i=0; i<this->NumParticle; i++)
		op.AddFrom(this->Velocity[i]);

	op.MultiplyFrom( 1.0/(double)(this->NumParticle) );

	for(size_t i=0; i<this->NumParticle; i++)
		this->Velocity[i].MinusFrom(op);
}




#ifdef DEBUG
const double power1=0.04;
const double coeff1=20;
#else
const double power1=0.04;
const double coeff1=0.2;
#endif


bool ParticleMDAnnealing::StopSignal=false;

////////////////////////////////////////////////////////////////////////////////////////
//ParticleMDAnnealing::ParticleMDAnnealing(DimensionType Dimension, size_t NumParticle, double InitialDensity, int RandomSeed, double Pressure, Potential & pot, std::ostream * output, double CoolingScheduleRescale, double MinDistance, bool UseSortedList) : Output(output), gen(RandomSeed), MinDistance(MinDistance)
//{
//	this->alpha=1.0- coeff1*std::pow( power1, Dimension)/CoolingScheduleRescale;
//	if(alpha==1)
//	{
//		std::cerr<<"Error in ParticleSimulatedAnnealling::ParticleSimulatedAnnealing : System is too complex, unable to find a cooling schedule!";
//		assert(false);
//	}
//	this->Tcut=0.00000100*std::pow(0.1, Dimension);
//
//	this->Tinit=100000;
//	//this->CellMoveProbability=1/(5.0*NumParticle-4.0);
//	this->CellMoveProbability=0;
//
//	this->Pressure=Pressure;
//	size_t ExpectedStep=std::log(Tinit/Tcut)/std::log(alpha);
//
//	//////////////////////
//	//initialize this->pList
//	if(NumParticle==0)
//	{
//		std::cerr<<"Error in MD_System::MD_System : NumParticle is 0\n";
//		assert(false);
//	}
//	unsigned short cellrank = std::floor( std::pow(NumParticle/ParticlePerCell, static_cast<double>(1)/Dimension) );//about 1 particles per cell
//	double initsize=std::pow(NumParticle/InitialDensity, static_cast<double>(1)/Dimension);
//	std::vector<GeometryVector> base;
//	for(DimensionType i=0; i<Dimension; i++)
//	{
//		base.push_back(GeometryVector(Dimension));
//		base.back().x[i]=initsize;
//		for(DimensionType j=i+1; j<Dimension; j++)
//			base.back().x[j]=initsize*(gen.RandomDouble()-0.5);
//	}
//	this->pList=new Configuration(Dimension, &base[0], initsize/cellrank, UseSortedList);
//	size_t trialcount=0;
//	for(size_t i=0; i<NumParticle;)
//	{
//		GeometryVector temp;
//		for(DimensionType j=0; j<Dimension; j++)
//			temp.x[j]=gen.RandomDouble();
//		bool TooNearNeighbor=false;
//		this->pList->IterateThroughNeighbors(temp, MinDistance, [&TooNearNeighbor](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
//		{
//			TooNearNeighbor=true;
//		}, &TooNearNeighbor);
//		if(TooNearNeighbor==false)
//		{
//			this->pList->Insert("C", temp);
//			i++;
//		}
//		trialcount++;
//		if(trialcount%100000==99999)
//		{
//			GeometryVector tbase[ ::MaxDimension];
//			for(DimensionType j=0; j<Dimension; j++)
//				tbase[j]=this->pList->GetBasisVector(j)*2;
//			this->pList->ChangeBasisVector(tbase);
//		}
//	}
//
//	pot.SetConfiguration(*this->pList);
//	double Enthalpy=pot.Energy()+Pressure*this->pList->PeriodicVolume();
//
//	///////////////////////
//	//initialize this->pSys
//	double TimeStep;
//	this->pSys = new ParticleMolecularDynamics(this->pList, TimeStep, 1.0);
//
//	if(this->Output!=nullptr)
//	{
//		(*this->Output).precision(14);
//		(*this->Output)<<"Simulated Annealing Started, alpha="<<alpha<<", T_init="<<Tinit<<", T_cut="<<Tcut<<", Cell Move Probability="<<CellMoveProbability<<", Pressure="<<Pressure<<", Mindistance="<<this->MinDistance<<'\n';
//	}
//	this->RelaxStructureFunction= [] (Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance) ->void
//	{
//		RelaxStructure_GSL(List, pot, Pressure, Switch, MinDistance, 10000);
//	};
//}
ParticleMDAnnealing::ParticleMDAnnealing(const Configuration & InitConfig, int RandomSeed, double Pressure, Potential & pot, bool BoxMCMove, std::ostream * output, double CoolingScheduleRescale, bool UseSortedList) : 
	Output(output), gen(RandomSeed), BoxMCMove(BoxMCMove), AutoTimeStep(true)
{
	if(BoxMCMove)
		this->RelaxStructureFunction= [] (Configuration & List, Potential & pot, double Pressure, double MinDistance) ->void
	{
		RelaxStructure_ConjugateGradient(List, pot, Pressure, 0, MinDistance, 10000);
	};
	else
		this->RelaxStructureFunction= [] (Configuration & List, Potential & pot, double Pressure, double MinDistance) ->void
	{
		RelaxStructure_ConjugateGradient(List, pot, Pressure, 0, MinDistance, 10000);
	};
	this->RelaxStructureBeforeAnnealing = true;

	this->SameLocalMinimaDetection=true;

	this->alpha=1.0- coeff1*std::pow( power1, InitConfig.GetDimension())/CoolingScheduleRescale;
	if(alpha==1)
	{
		std::cerr<<"Error in ParticleSimulatedAnnealling::ParticleSimulatedAnnealing : System is too complex, unable to find a cooling schedule!";
		assert(false);
	}
	this->Tcut=0.00000100*std::pow(0.1, InitConfig.GetDimension());

	this->Tinit=1000;
	//this->CellMoveProbability=1/(5.0*InitConfig.NumParticle()-4.0);

	this->Pressure=Pressure;
	size_t ExpectedStep=std::log(Tinit/Tcut)/std::log(alpha);


	///////////////////////
	//initialize this->pSys
	double LengthScale = std::pow(InitConfig.PeriodicVolume() / InitConfig.NumParticle(), 1.0 / (InitConfig.GetDimension()));
	double defaultTimeStep = 0.01*LengthScale / std::sqrt(Tinit / 2.0);
	this->pSys = new ParticleMolecularDynamics(InitConfig, defaultTimeStep, 1.0);
	//////////////////////
	//initialize this->pList
	this->pList=&this->pSys->Position;
	pot.SetConfiguration(*this->pList);
	pot.Energy();
	//double Enthalpy=pot.Energy()+Pressure*this->pList->PeriodicVolume();

	this->Continue = false;
	this->ConfigPackName = "MDAnneal";
}


#include "MonteCarlo.h"
class StopMDAnnealing
{
public:
	StopMDAnnealing()
	{}
};

const double ZeroMoveStopTolerence = 1e-10;
void ParticleMDAnnealing::Anneal(Potential & pot, std::function<void (void)> Callback)
{
	if(this->Output!=nullptr)
	{
		//(*this->Output).precision(14);
		(*this->Output)<<"MD Annealing Started, alpha="<<alpha<<", T_init="<<Tinit<<", T_cut="<<Tcut<<", Pressure="<<Pressure<<'\n';
	}

	//this->pSys->Evolve(10000, pot);
	this->pSys->SetRandomSpeed(this->Tinit, gen);

	double p = this->Pressure;
	std::ostream * op= Output;
	double prevT=-1.0;

	ConfigurationPack OutputPack(ConfigPackName);
	if (Continue && OutputPack.NumConfig()>0)
		(*this->pList) = OutputPack.GetConfig(OutputPack.NumConfig() - 1);
	else
	{
		OutputPack.Clear();
		if(RelaxStructureBeforeAnnealing)
			RelaxStructureFunction(*this->pList, pot, Pressure, 0.0);
	}

	if(this->pSys != nullptr)
	{
		pot.SetConfiguration(*this->pList);
		//this->pSys->RecordEnergy=pot.Energy()+p*this->pList->PeriodicVolume();

		//ThermoCool cool2(this->alpha, this->Tinit, this->Tcut, 0, 20, 5000, 100);
		ExpCool cool2(this->alpha, this->Tinit, this->Tcut, 0);
		double prevCheckRelaxStructureT=this->Tinit;
		double prevCheckRelaxStructureH=pot.Energy()+p*this->pList->PeriodicVolume();
		size_t GotSameLocalCount=0;
		bool * pStopSignal=& ParticleMDAnnealing::StopSignal;

		CellMove * pCellMove = nullptr;
		if(BoxMCMove)
			pCellMove = new CellMove(Pressure);

		try
		{
			while(cool2.Continue==true)
			{
				//debug temp 
				//unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

				double Temperature = cool2.Temperature;

				for(double i=0; i<cool2.NumTry; i+=100000*this->pSys->TimeStep)
				{
					this->pSys->AndersonEvolve(20, pot, Temperature, 0.1, gen);  
					if(BoxMCMove)
					{
						pot.SetConfiguration(this->pSys->Position);
						double CurrentH=pot.Energy()+p*this->pList->PeriodicVolume();
						for(int i=0; i<10; i++)
						{
							double PreFactor=1.0;
							double dE=pCellMove->DeltaEnergy(this->pSys->Position, &pot, this->gen, PreFactor, false, CurrentH);
							if(gen.RandomDouble() < PreFactor*std::exp((-1)*dE/Temperature))
							{
								pCellMove->Accept(this->pSys->Position, &pot);
								CurrentH+=dE;
							}
						}
					}
				}

				//::Output("MDAnneal_CurrentStructure", *this->pList);

				pot.SetConfiguration(*this->pList);
				double Enthalpy=pot.Energy()+p*this->pList->PeriodicVolume();

				//debug temp 
				//fp_control_state = _controlfp(_MCW_EM, _MCW_EM);

				cool2.Report(Enthalpy);

				if (op != nullptr && (Temperature != prevT || Verbosity >4))
				{
					(*op) << std::time(nullptr) - ::ProgramStart << " \t" << this->pSys->Time << " \t" << this->pSys->TimeStep << " \t" << Temperature << " \t" << this->pList->PeriodicVolume() << " \t" << Enthalpy << " \t" << this->pSys->GetKineticEnergy() << " \t" << OutputPack.NumConfig() - 1 << '\n';
					op->flush();
					OutputPack.AddConfig(*this->pList);
				}

				if (Temperature != prevT)
				{
					Callback();

					prevT=Temperature;

					if(this->AutoTimeStep)
						//this->pSys->Evolve_AutoTimeStep(50, pot, 5e-6);//soft potentials (stealthy, targetS(0), ...)
						this->pSys->Evolve_AutoTimeStep(50, pot, 1e-7);//soft potentials (stealthy, targetS(0), ...), conservative
						//this->pSys->Evolve_AutoTimeStep(50, pot, 1e-7);//stiff-core potentials (Lennard-Jones, R^n, ...), not always working
						//this->pSys->Evolve_AutoTimeStep(50, pot, 1e-8);//stiff-core potentials (Lennard-Jones, R^n, ...), more conservative
				}
				if((*(pStopSignal))==true)
					throw StopMDAnnealing();

				if( (Temperature<0.5*prevCheckRelaxStructureT || (Temperature<0.85*prevCheckRelaxStructureT&&GotSameLocalCount>0))&& this->SameLocalMinimaDetection )
				{
					//check if we get the same enthalpy after Relaxing
					Configuration temp(*this->pList);

					this->RelaxStructureFunction(temp, pot, p, 0.0);

					//debug temp
					//::PlotDirectionalStructureFactor2D(temp, 30, "CurrentInherent_Sk");

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
							throw StopMDAnnealing();
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
				pot.SetConfiguration(*this->pList);
			}
		}
		catch(StopMDAnnealing & a)
		{}
		if(pCellMove!=nullptr)
			delete pCellMove;
	}

	if(ParticleMDAnnealing::StopSignal==false)
	{
		pot.SetConfiguration(*this->pList);
		if(op!=nullptr)
		{
			(*op)<<"Optimization Start, Starting size:"<<this->pList->PeriodicVolume()<<", Enthalpy:"<<pot.Energy()+this->Pressure*this->pList->PeriodicVolume()<<'\n';
		}
		this->RelaxStructureFunction(*this->pList, pot, p, 0.0);
		pot.SetConfiguration(*this->pList);
		if(op!=nullptr)
		{
			(*op)<<"Optimization Finish, Ending size:"<<this->pList->PeriodicVolume()<<", Enthalpy:"<<pot.Energy()+this->Pressure*this->pList->PeriodicVolume()<<'\n';
		}
	}
}


