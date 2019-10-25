#include "Optimize.h"
#include "LatticeSumSolver.h"
#include "withGaussianInvParameterization.h"
#include "Plots.h"
#include "LJPotential.h"
#include "MC_System.h"
#include "PhononFrequency.h"
#include "Polynomial.h"
#include "Interfaces.h"
#include "PairCorrelation.h"
#include "StructureFactor.h"
#include "Degeneracy.h"
#include "ThreeBody.h"
#include "StructureOptimization.h"
#include "Solvers.h"
#include "MD_System.h"

const size_t LargeMCNumParticles=24;
const size_t LargeMCNumConfigs=4;

const size_t FixBoxMCSystemSideLengthFactor=2;
const size_t FixBoxMCNumConfigs=4;

#include <iostream>
#include <sstream>
void RelaxVacancy(Configuration list, Potential & pot, size_t MultiplicateNumber, double Pressure)
{
	if(MultiplicateNumber%2==0)
		std::cerr<<"Warning in Interfaces.cpp : RelaxVacancy : MultiplicateNumber is odd, the deleted particle may not be in the center!\n";
	Configuration list2=MultiplicateStructure(list, MultiplicateNumber);
	RelaxStructure(list2, pot, Pressure, 0.0);
	std::fstream ofile("RelaxVacancy_output.txt", std::fstream::out);
	ofile.precision(14);

	pot.SetConfiguration(list2);
	double OriginalH=pot.Energy()+list2.PeriodicVolume()*Pressure;
	ofile<<"original H/N:"<<(OriginalH)/list2.NumParticle()<<'\n';

	double OriginalHofNMinusOneParticles=OriginalH/list2.NumParticle()*(list2.NumParticle()-1);

	//delete the particle in the middle to make plots look better
	//if the Multiplicate time is odd, then this is the middle atom
	//todo : add "find nearest" feature in Configuration so that even if Multiplicate time is even, 
	//we can still find the center particle to delete
	//Configuration::particle * todelete=list2.GetParticle(list2.NumParticle()/2);
	//list2.DeleteParticle(todelete);
	list2.DeleteParticle(list2.NumParticle()/2);

	Output("RelaxVacancy_Before", list2);
	Plot("RelaxVacancy_Before", list2);
	pot.SetConfiguration(list2);
	ofile<<"before relaxation H increased:"<<(pot.Energy()+list2.PeriodicVolume()*Pressure)-OriginalHofNMinusOneParticles<<'\n';

	RelaxStructure(list2, pot, Pressure, 0.0);

	Output("RelaxVacancy_after", list2);
	Plot("RelaxVacancy_after", list2);
	pot.SetConfiguration(list2);
	ofile<<"after relaxation H increased:"<<(pot.Energy()+list2.PeriodicVolume()*Pressure)-OriginalHofNMinusOneParticles<<'\n';
}

void OutputThetaSeries(LatticeSumSolver & solver, std::ostream & output)
{
	if(solver.Terms.size()<2)
	{
		output<<"Warning in OutputThetaSeries : No theta series!\n";
		return;
	}
	output<<"Theta Series:\n";
	for(auto iter=solver.Terms.begin()+1; iter!=solver.Terms.end(); iter++)
	{
		output<<(double)(iter->number)/solver.Terms[0].number<<" \t"<<iter->distance*iter->distance/solver.Terms[1].distance/solver.Terms[1].distance<<'\n';
	}
	output<<"Neighbor Series:\n";
	for(auto iter=solver.Terms.begin(); iter!=solver.Terms.end(); iter++)
	{
		output<<iter->number<<" \t"<<iter->distance/*/solver.Terms[1].distance*/<<'\n';
	}
}
void OutputCumulativeCoordination(LatticeSumSolver & solver, std::ostream & output)
{
	if(solver.Terms.size()<2)
	{
		output<<"Warning in OutputThetaSeries : No theta series!\n";
		return;
	}
	output<<"r\t Z(r)\n";
	double sum=0.0;
	for(auto iter=solver.Terms.begin()+1; iter!=solver.Terms.end(); iter++)
	{
		output<<(iter->distance/solver.Terms[1].distance)*(1-1e-8)<<" \t"<<sum<<'\n';
		sum+=(double)(iter->number)/solver.Terms[0].number;
		output<<iter->distance/solver.Terms[1].distance<<" \t"<<sum<<'\n';
	}
}
void OutputCoordinationNumber(LatticeSumSolver & solver, std::ostream & output)
{
	if(solver.Terms.size()<2)
	{
		output<<"Warning in OutputThetaSeries : No theta series!\n";
		return;
	}
	output<<"r \t Z(r):\n";
	for(auto iter=solver.Terms.begin()+1; iter!=solver.Terms.end(); iter++)
	{
		output<<iter->distance/solver.Terms[1].distance<<" & "<<(double)(iter->number)/solver.Terms[0].number<<"\n";
	}
}

long RedundentSolver(std::vector<LatticeSumSolver *> solvers, LatticeSumSolver & tested)
{
	for(auto iter=solvers.begin(); iter!=solvers.end(); iter++)
	{
		bool result=SameSolver(*iter, &tested);
		if(result == true)
			return iter-solvers.begin();
	}
	return -1;
}

void RestoreStatus(ParameterSet & param, std::vector<LatticeSumSolver * > & solver, size_t & NumParam, size_t & NumCompetitor, std::string Prefix = "")
{
	for(NumParam =0; ; NumParam++)
	{
		std::stringstream stemp;
		char PrefixBuffer[64];
		stemp<<Prefix<<"Optimization_"<<NumParam<<"_Parameters.txt";
		stemp.getline(PrefixBuffer, 63);
		std::fstream infile(PrefixBuffer, std::fstream::in);
		if(infile.good()==false)
			break;
		param.Read(infile);
	}
	for(NumCompetitor =0; ; NumCompetitor++)
	{
		std::stringstream stemp;
		char PrefixBuffer[64];
		stemp<<Prefix<<"Annealing_"<<NumCompetitor<<".pos";
		stemp.getline(PrefixBuffer, 63);
		std::fstream infile(PrefixBuffer, std::fstream::in);
		if(infile.good()==false)
			break;
		infile.close();
		solver.push_back(new LoadStructure(PrefixBuffer));
	}
}


#include <omp.h>
int InvThermoIteration(std::vector<LatticeSumSolver *> & solver, ParameterSet * & pParam, const size_t ReOptimizeTime, const std::vector<size_t> & NumInstances, size_t AnnealTargetNumRequirement, void (*PlotPhononFunction)(const Configuration &, PairPotential &, char *, double, size_t), double CoolSchedule = 1.0, double * PressureOutput=nullptr, int seed=0, bool AllowSkipOptimization=true)
{
	::Output("Target", solver[0]->GetStructure());
	std::fstream TargetTheta("Target.ThetaSeries", std::fstream::out);
	::OutputThetaSeries(*solver[0], TargetTheta);
	TargetTheta.close();
	size_t NumPosOutput=0;
	size_t NumParam=0;
	size_t TargetAnnealed=0;
	RandomGenerator gen(seed);
	size_t NumMCThreads = ::MCNumParallelConfigurations ;
	if(NumMCThreads == 0)
		NumMCThreads = omp_get_num_procs();

	assert(solver.size()>0);
	const Configuration Target=solver[0]->GetStructure();
	const double TargetStructureSpecificVolume=Target.PeriodicVolume()/Target.NumParticle();
	::DimensionType dimension = Target.GetDimension();
	RestoreStatus(*pParam, solver, NumParam, NumPosOutput);
	std::cout<<"Loaded "<<NumParam<<"th parameter, and "<<NumPosOutput<<"Competitors\n";
	logfile<<"Loaded "<<NumParam<<"th parameter, and "<<NumPosOutput<<"Competitors\n";
	const size_t MaxNumParticle=NumInstances.size(); 
	volatile bool CompetitorFound = false;
	volatile bool ContinueSimulatedAnnealing = true;
	for(;;)
	{
		std::cout<<"Original Parameterization:";
		pParam->Write(std::cout);
		double Result=-1;
		double CenterPressure=0;
		double DeltaPressure=0;

		IsotropicPairPotential * pPotential=pParam->GetPotential();
		CenterPressure=GetPressure(solver[0], TargetStructureSpecificVolume, *pPotential);
		bool TargetLocallyStable=EnthalpyDifference::LocallyStable(Target, *pPotential, CenterPressure, pParam->MinDistance());
		//bool TargetLocallyStable=true;
		delete pPotential;

		if(solver.size()<2)
		{
			std::cout<<"No known Competitor, skip optimization!\n";
			logfile<<"No known Competitor, skip optimization!\n";
		}
		else
		{
			double HDiff=EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, CenterPressure);
			if(HDiff<0.0 && TargetLocallyStable && AllowSkipOptimization)
			{
				std::cout<<"Target is stable, skip optimization!\n";
				logfile<<"Target is stable, skip optimization!\n";
			}
			else
			{

				std::cout<<"Start Optimization.\n";
				logfile<<"Start Optimization.\n";
				bool OptimizeResult=ReOptimize(pParam, solver, ReOptimizeTime, TargetStructureSpecificVolume, Result, NumParam, CenterPressure, DeltaPressure, gen);

				LowestFrequencySquared = 0.0;
				if(PlotPhononFunction != nullptr)
				{
					std::stringstream stemp;
					char PrefixBuffer[32];
					stemp<<"Optimization_"<<NumParam-1<<"_Phonons";
					stemp.getline(PrefixBuffer, 31);
					IsotropicPairPotential * p=pParam->GetPotential();
					Configuration s=solver[0]->GetStructure();
					PlotPhononFunction(s, *p, PrefixBuffer, 1.0, 100);
					delete p;
				}
				if(OptimizeResult==false)
				{
					std::cout<<"Can't find a potential! Exiting\n";
					logfile<<"Can't find a potential! Exiting\n";
					break;
				}
				std::cout<<'\n';
				logfile<<'\n';
			}
		}

		TargetAnnealed = 0;
		double LPressure=CenterPressure-DeltaPressure;
		double HPressure=CenterPressure+DeltaPressure;
		double HDiff;
		
		if(solver.size()>=2)
		{
			HDiff=EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, CenterPressure);
			//double checks
			if(EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, HPressure)>0.0)
			{
				std::cerr<<EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, HPressure);
				assert(false);
			}
			if(EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, LPressure)>0.0)
			{
				std::cerr<<EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, LPressure);
				assert(false);
			}
		}

		bool HasImaginaryPhononFrequency = false;
		if(LowestFrequencySquared<-0.1)
		{
			LowestFrequencySquared=0;
			HasImaginaryPhononFrequency=true;
		}

		double pMax=(-1)*Result;
		std::vector<ParticleSimulatedAnnealing *> vpSimuAnneal;
anneal:
		IsotropicPairPotential * pPot = (*pParam).GetPotential();
		Configuration Obj=solver[0]->GetStructure();
		if(solver.size()>=2)
		{
			std::cout<<"dH:"<<HDiff<<", Pressure:"<<LPressure<<", Lowest Eignevalue of Deformation Tensor:"<<LowestEigenValueOfDeformations(Obj, *pPot, CenterPressure)<<", Lowest Eignevalue of Stiffness Tensor:"<<LowestEigenValueOfStiffness(Obj, *pPot, CenterPressure)<<'\n';
			logfile<<"dH:"<<HDiff<<", Pressure:"<<LPressure<<", Lowest Eignevalue of Deformation Tensor:"<<LowestEigenValueOfDeformations(Obj, *pPot, CenterPressure)<<", Lowest Eignevalue of Stiffness Tensor:"<<LowestEigenValueOfStiffness(Obj, *pPot, CenterPressure)<<'\n';
		}
		//start simulated annealing
		std::cout<<"Start Simulated Annealing\n";
		logfile<<"Start Simulated Annealing\n";

		size_t ObjN=solver[0]->GetStructure().NumParticle();


		std::vector<double> Pressures;
		Pressures.push_back(CenterPressure);
		for(auto iter=Pressures.begin(); iter!=Pressures.end(); iter++)
		{
			double Vol=TargetStructureSpecificVolume;
			double Ene = solver[0]->LatticeSum( Vol, *pPot);
			double Ent = Ene + Vol*(*iter);
			std::cout<<"at P="<<*iter<<", Target V/N="<<Vol<<", E/N="<<Ene<<", H/N="<<Ent<<'\n';
			logfile<<"at P="<<*iter<<", Target V/N="<<Vol<<", E/N="<<Ene<<", H/N="<<Ent<<'\n';
		}

		vpSimuAnneal.clear();
		if(PressureOutput!=nullptr)
			(*PressureOutput)=CenterPressure;

		for(size_t i=1; i<= MaxNumParticle ; i++)
		{
			for(size_t j=0; j<NumInstances[i-1]; j++)
			{
				for(auto iter=Pressures.begin(); iter!=Pressures.end(); iter++)
					vpSimuAnneal.push_back(new ParticleSimulatedAnnealing(dimension, i, 1, gen.RandomDouble()*65535, *iter, *pPot, nullptr, CoolSchedule, pParam->MinDistance(), 1));
			}
		}
		if(HasImaginaryPhononFrequency)
		{
			std::cout<<"Imaginary Phonons detected, add large configurations in simulated annealing\n";
			logfile<<"Imaginary Phonons detected, add large configurations in simulated annealing\n";
			for(size_t i=0; i<LargeMCNumConfigs; i++)
				vpSimuAnneal.push_back(new ParticleSimulatedAnnealing(dimension, LargeMCNumParticles, 1, gen.RandomDouble()*65535, CenterPressure, *pPot, nullptr, CoolSchedule*0.3, pParam->MinDistance(), 1));
		}
		for(size_t i=0; i<FixBoxMCNumConfigs; i++)
		{
			Configuration temp=MultiplicateStructure(solver[0]->GetStructure(), FixBoxMCSystemSideLengthFactor);
			vpSimuAnneal.push_back(new ParticleSimulatedAnnealing(dimension, temp.NumParticle(), 1, gen.RandomDouble()*65535, CenterPressure, *pPot, &logfile, CoolSchedule*2.0, pParam->MinDistance(), 1, true));
			vpSimuAnneal.back()->CellMoveProbability=0.0;
			(*vpSimuAnneal.back()->pList) = temp;
		}

		CompetitorFound = false;
		ContinueSimulatedAnnealing = true;
		ParticleSimulatedAnnealing::StopSignal=false;
		signed long end=vpSimuAnneal.size();


		std::vector<LatticeSumSolver *> Recents;

#pragma omp parallel default(shared) num_threads(NumMCThreads)
		{
			IsotropicPairPotential * ThreadpPot = (*pParam).GetPotential();
#pragma omp for schedule(dynamic)
			for(signed long i=0; i<end; i++)
			{
				if(CompetitorFound==true)
				{
					delete vpSimuAnneal[i];
					continue;
				}
				vpSimuAnneal[i]->Anneal(*ThreadpPot);
#pragma omp critical
				{
					size_t N=vpSimuAnneal[i]->pList->NumParticle();
					double VoverN=vpSimuAnneal[i]->pList->PeriodicVolume()/N;
					double EoverN=ThreadpPot->Energy(*vpSimuAnneal[i]->pList)/N;
					double HoverN=EoverN+vpSimuAnneal[i]->Pressure*VoverN;
					Output(std::string("temp"),*vpSimuAnneal[i]->pList);
					LoadStructure tempSolver("temp.pos");
					std::cout<<"at time"<<std::time(nullptr)-ProgramStart<<", "<<i<<"th, N="<<N<<", P="<<vpSimuAnneal[i]->Pressure<<", V/N="<<VoverN<<", E/N="<<EoverN<<", H/N="<<HoverN<<", D="<<SolverDistance(solver[0], &tempSolver)<<", ";
					logfile<<"at time"<<std::time(nullptr)-ProgramStart<<", "<<i<<"th, N="<<N<<", P="<<vpSimuAnneal[i]->Pressure<<", V/N="<<VoverN<<", E/N="<<EoverN<<", H/N="<<HoverN<<", D="<<SolverDistance(solver[0], &tempSolver)<<", ";
					bool Threatening=false;

					//test the annealing result
					volatile long stemp=RedundentSolver(solver, tempSolver);
					volatile long stemp2=RedundentSolver(Recents, tempSolver);
					if(stemp!=(-1))
					{
						std::cout<<"same as:"<<solver[stemp]->Tag()<<", ";
						logfile<<"same as:"<<solver[stemp]->Tag()<<", ";
						if(stemp==0)
							TargetAnnealed++;
					}
					else if(stemp2!=(-1))
					{
						std::cout<<"same as:"<<Recents[stemp2]->Tag()<<", ";
						logfile<<"same as:"<<Recents[stemp2]->Tag()<<", ";
					}
					//else
					{
						solver.push_back(&tempSolver);
						if(EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, vpSimuAnneal[i]->Pressure)>=0)
						{
							std::cout<<"at P="<<vpSimuAnneal[i]->Pressure<<", dH="<<EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, vpSimuAnneal[i]->Pressure)<<", thus threatening.";
							logfile<<"at P="<<vpSimuAnneal[i]->Pressure<<", dH="<<EnthalpyDifference::EnthalpyDifference(*pParam, solver, TargetStructureSpecificVolume, vpSimuAnneal[i]->Pressure)<<", thus threatening.";
							Threatening=true;
						}
						solver.pop_back();
						if(Threatening==false)
						{
							std::cout<<"not threatening\n";
							logfile<<"not threatening\n";
						}
					}

					if(Threatening /*&& (stemp==(-1)) && (stemp2==(-1)) */&& ParticleSimulatedAnnealing::StopSignal==false)
					{
						//threatening Competitor found
						CompetitorFound = true;
						ContinueSimulatedAnnealing = false;
						ParticleSimulatedAnnealing::StopSignal=true;
						std::cout<<"Write as Annealing_"<<NumPosOutput<<"...";
						logfile<<"Write as Annealing_"<<NumPosOutput<<"...";
						std::cout.flush();

						{
							std::stringstream stemp;
							char PrefixBuffer[32];
							stemp<<"Annealing_"<<NumPosOutput;
							stemp.getline(PrefixBuffer, 31);
							std::string temp(PrefixBuffer);
							Output(temp, *vpSimuAnneal[i]->pList);
							Plot(temp, *vpSimuAnneal[i]->pList);
							temp=temp+".pos";
							solver.push_back(new LoadStructure(temp.c_str()));
						}

						{
							std::stringstream stemp;
							char PrefixBuffer[32];
							stemp<<"Annealing_"<<NumPosOutput;
							stemp.getline(PrefixBuffer, 31);
							::Plots((*pParam), & solver, PrefixBuffer, solver.size()-1);
						}

						{
							std::stringstream stemp;
							char PrefixBuffer[32];
							stemp<<"Annealing_"<<NumPosOutput<<".ThetaSeries";
							stemp.getline(PrefixBuffer, 31);
							std::fstream outfile2(PrefixBuffer, std::fstream::out);
							OutputThetaSeries(*solver.back(), outfile2);
						}

						std::cout<<"done\n";
						logfile<<"done\n";
						NumPosOutput++;

						Recents.push_back(solver.back());
						solver.pop_back();
					}
					delete vpSimuAnneal[i];
				}
			}
			delete ThreadpPot;
		}
		std::cout<<'\n';
		logfile<<'\n';

		for(auto iter=Recents.begin(); iter!=Recents.end(); iter++)
			solver.push_back(*iter);
		Recents.clear();
		delete pPot;


		if(ContinueSimulatedAnnealing)
		{
			std::cout<<"No threatening Competitor found in simulated annealing, PERMANANTLY INCREASE COOLING SCHEDULE LENGTH!\n";
			logfile<<"No threatening Competitor found in simulated annealing, PERMANANTLY INCREASE COOLING SCHEDULE LENGTH!\n";
			CoolSchedule*=2;
			if(TargetAnnealed>AnnealTargetNumRequirement)
				break;
			else
				goto anneal;
		}
	}
	return 0;
}
int CheckResult(std::vector<LatticeSumSolver *> & solver, ParameterSet * & pParam, const size_t MCNumParticle, const size_t MCNumConfig, void (*PlotPhononFunction)(const Configuration &, PairPotential &, char *, double, size_t), double CoolSchedule = 1.0)
{
	::AlsoWriteEPS=true;
	size_t temp1, temp2;
	RestoreStatus(*pParam, solver, temp1, temp2);
	std::cout<<"Verifying results!\nParameters are:";
	logfile<<"Verifying results!\nParameters are:";
	pParam->Write(std::cout);
	pParam->Write(logfile);
	std::fstream paramfile("CheckResult_Parameters.txt", std::fstream::out);
	pParam->Write(paramfile);
	std::cout<<'\n';
	logfile<<'\n';
	IsotropicPairPotential * pPot = pParam->GetPotential();

	::Plots(pPot, &solver, "CheckResult");

	auto stru=solver[0]->GetStructure();
	double MCPressure=GetPressure(solver[0], stru.PeriodicVolume()/stru.NumParticle(), *pPot);
	PlotPhononFunction(stru, *pPot, "CheckResult", 1.0, 500);

	RelaxVacancy(stru, *pPot, 11, MCPressure);

	//print elastic constants
	std::fstream ofile2("CheckResult_Target.ElasticConstants", std::fstream::out);
	PrintAllElasticConstants(ofile2, stru, *pPot, MCPressure);
	ofile2.close();

	std::cout<<"Do Simulated Annealing for the final Verification!\n";
	logfile<<"Do Simulated Annealing for the final Verification!\n";
	std::cout<<"Num of Particle="<<MCNumParticle<<", Pressure="<<MCPressure<<". For "<<MCNumConfig<<" Configurations.\n";
	logfile<<"Num of Particle="<<MCNumParticle<<", Pressure="<<MCPressure<<". For "<<MCNumConfig<<" Configurations.\n";
	::DimensionType dim=solver[0]->GetStructure().GetDimension();
	signed long NumConfig=MCNumConfig;

#pragma omp parallel for
	for(signed int i=0; i<NumConfig; i++)
	{
		IsotropicPairPotential * ThreadpPot = pParam->GetPotential();
		ParticleSimulatedAnnealing a(dim, MCNumParticle, 1.0, i, MCPressure, *ThreadpPot, i==0?&logfile:nullptr, CoolSchedule, pParam->MinDistance());
		a.Anneal(*ThreadpPot);
		std::stringstream Number;
		Number<<i;
		std::string name("CheckResult_MCAnneal_");
		std::string name2("");
		Number>>name2;
		name+=name2;
		Output(name.c_str(), *a.pList);
		std::cout<<i<<"th result completed! E/N="<<ThreadpPot->Energy(*a.pList)/MCNumParticle<<", V/N="<<a.pList->PeriodicVolume()/MCNumParticle<<", H/N="<<(ThreadpPot->Energy(*a.pList)+MCPressure*a.pList->PeriodicVolume())/MCNumParticle<<'\n';
		logfile<<i<<"th result completed! E/N="<<ThreadpPot->Energy(*a.pList)/MCNumParticle<<", V/N="<<a.pList->PeriodicVolume()/MCNumParticle<<", H/N="<<(ThreadpPot->Energy(*a.pList)+MCPressure*a.pList->PeriodicVolume())/MCNumParticle<<'\n';

		std::cout<<"Calculate Theta Series\n";
		logfile<<"Calculate Theta Series\n";
		LoadStructure l((name+std::string(".pos")).c_str());
		LJPotential lj(dim, 1, 1, 4);
		l.LatticeSum(1.0, lj);
		std::fstream ofile((name+std::string(".ThetaSeries")).c_str(), std::fstream::out);
		OutputThetaSeries(l, ofile);

		Plot(name, *a.pList);
		delete ThreadpPot;
	}

	delete pPot;
	::AlsoWriteEPS=false;

	return 0;
}

//the first version, all parameters hard-coded
namespace
{

	double Cutoff=16.1;
	double core=0.75;
	class ErdalPotential : public IsotropicPairPotential
	{
	public:
		double W0, V0, d, Sigma, K;
		ErdalPotential() : IsotropicPairPotential(2, Cutoff), W0(1.6), V0(1.0), d(1.0), Sigma(core), K(1.0)
		{
		}
		virtual double IsotropicEnergyFormula(double distance, const AtomInfo & a1, const AtomInfo & a2)
		{
			double result= (distance<Sigma)? ::MaxEnergy : 0.0;
			result+=V0*std::exp((-1)*K*distance)/((1)*K*distance);

			if(distance < d)
				result+=(-1)*W0*(1.0-3.0*distance/2.0/d+distance*distance*distance/2.0/d/d/d);

			return result;
		}
	};

	class ErdalParameters: public ParameterSet
	{
	public:
		ErdalParameters(LatticeSumSolver * TargetStructure)
		{
			this->Parameters=new double [1];

			this->Parameters[0]=0.0;
		}
		ErdalParameters(const ErdalParameters & source):ParameterSet(source)
		{
		}
		virtual IsotropicPairPotential * GetPotential(void) const
		{
			return new ErdalPotential();
		}
		virtual double Penalty(void) const//calculate the penalty from OptimizableParameters
		{
			return 1;
		}
		virtual void GetBound(double * LowerBound, double * UpperBound) const//get range of each parameters
		{
			LowerBound[0]=0.0;
			UpperBound[0]=0.0;
		}
		virtual bool Evolve(void)//evolve to a more complex Parameterization, return true if succed, return false if can't evolve any more
		{
			return false;
		}
		virtual bool Iterate(void)//for some non-continuous parameter, iterate through it. Intended to be used when it has never evolved
		{
			return false;
		}
		virtual ParameterSet * clone() const
		{
			return new ErdalParameters(*this);
		}
		virtual void Write(std::ostream & out) const
		{
			out<<'\n';
		}
		virtual void Read(std::istream & in)
		{
		}
		virtual double MinDistance(void) const//when r<MinDistance, the potential should have no interesting behavior, guaranteed to fix unchanged during the life of this object
		{
			return core;
		}
		virtual double MaxDistance(void) const//the Rcut of the potential should not exceed this, guaranteed to fix unchanged during the life of this object
		{
			return Cutoff;
		}
		virtual size_t NumParameters(void) const
		{
			return 1;
		}
		virtual void SignalGlobalOptimization(void)//a signal that global optimization will be carried out
		{
		}
		virtual void SignalLocalOptimization(void)//vice versa
		{
		}
		virtual double SharpestRelativeSize(void)
		{
			return 0.02;
		}
	};
};



int InvStatMechCLI()
{
	int InvStatMechDebug();
	try
	{
		char tempstring[1000];
		std::istream & ifile=std::cin;
		std::ostream & ofile=std::cout;
		std::vector<LatticeSumSolver *> solver(1, nullptr);
		ParameterSet * pParam = nullptr;
		size_t ReOptimizeTime = 100;
		std::vector<size_t> NumInstances;
		size_t AnnealTargetNumRequirement = 1;
		size_t RandomSeed =0;
		double Pressure = -1.0;
		double CoolingScheduleRescale=0.0;
		bool AllowSkipOptimization=true;
		void (*PlotPhononFunction)(const Configuration &, PairPotential &, char *, double, size_t) = nullptr;

		for(;;)
		{
			ifile>>tempstring;
			if(strcmp(tempstring, "Exit")==0)
			{
				for(auto iter=solver.begin(); iter!=solver.end(); iter++)
					if(*iter != nullptr)
						delete *iter;
				if(pParam != nullptr)
					delete pParam;
				return 0;
			}
			else if(strcmp(tempstring, "Target")==0)
			{
				if(solver[0]!=nullptr)
					delete solver[0];
				solver[0]=nullptr;

				ifile>>tempstring;
				solver[0] = new LoadStructure(tempstring);
			}
			else if(strcmp(tempstring, "Competitor")==0)
			{
				ifile>>tempstring;
				solver.push_back(new LoadStructure(tempstring));
			}
			else if(strcmp(tempstring, "Competetor")==0)//backward support for versions containing misspellings
			{
				ifile>>tempstring;
				solver.push_back(new LoadStructure(tempstring));
			}
			else if(strcmp(tempstring, "Debug")==0)
			{
				InvStatMechDebug();
			}
			else if(strcmp(tempstring, "Parameterization")==0)
			{
				if(solver[0] == nullptr)
				{
					std::cerr<<"Specify Target before specifing Parameterization!\n";
				}
				else
				{
					if(pParam != nullptr)
						delete pParam;
					pParam = nullptr;

					ifile>>tempstring;
					if(strcmp(tempstring, "withGaussianInv")==0)
						pParam = new withGaussianInvParameters(solver[0]);
					if(strcmp(tempstring, "Erdal")==0)
						pParam = new ErdalParameters(solver[0]);
					if(strcmp(tempstring, "TunnelwithGaussianInv")==0)
					{
						//the parameterization specifically for tunnel FCC crystal
						pParam = new TunnelwithGaussianInvParameters();
						pParam->Evolve();
					}
					else if(strcmp(tempstring, "Polynomial")==0)
						pParam = new PolynomialParameterization(solver[0]);
					else if(strcmp(tempstring, "GaussianCore")==0)
						pParam = new GaussianCoreParameterization(solver[0]);
					//else if(strcmp(tempstring, "SoftCore")==0)
					//	pParam = new SoftCoreParameterization(solver[0]);
					else if(strcmp(tempstring, "Combined")==0)
					{
						ParameterSet * p1=new withGaussianInvParameters(solver[0]);
						p1->Evolve();
						ParameterSet * p2=new GaussianCoreParameterization(solver[0]);
						pParam = new CombinedParameterization(p1, p2);
					}
					else
						std::cerr<<"Unrecongnized Parameterization!\n";
				}
			}
			else if(strcmp(tempstring, "ReOptimizeTime")==0)
			{
				ifile>>ReOptimizeTime;
			}
			else if(strcmp(tempstring, "RandomSeed")==0)
			{
				ifile>>RandomSeed;
			}
			else if(strcmp(tempstring, "MCNumInstances")==0)
			{
				while(ifile.peek()!='\n')
				{
					size_t temp;
					ifile>>temp;
					if(ifile.fail() == true)
					{
						ifile.clear();
						break;
					}
					NumInstances.push_back(temp);
				}
			}
			else if(strcmp(tempstring, "AchieveAtLeast")==0)
			{
				ifile>>AnnealTargetNumRequirement;
			}
			else if(strcmp(tempstring, "MCNumParallelConfigurations")==0)
			{
				ifile>> ::MCNumParallelConfigurations;
			}
			else if(strcmp(tempstring, "CoolingSchedule")==0)
			{
				ifile>> CoolingScheduleRescale;
			}
			else if(strcmp(tempstring, "Symmetry")==0)
			{
				if(solver[0] == nullptr)
				{
					std::cerr<<"Specify Target before specifing symmetry!\n";
				}
				else
				{
					ifile>>tempstring;
					if(strcmp(tempstring, "Hexagonal")==0)
					{
						if(solver[0]->GetStructure().GetDimension()==2)
							PlotPhononFunction = & PlotPhononFrequencies_2DHexagonal;
						else if(solver[0]->GetStructure().GetDimension()==3)
							PlotPhononFunction = & PlotPhononFrequencies_3DHexagonal;
						else
							std::cerr<<"Unrecongnized symmetry!\n";
					}
					else if(strcmp(tempstring, "Cubic")==0 || strcmp(tempstring, "Square")==0 || strcmp(tempstring, "HyperCubic")==0)
					{
						if(solver[0]->GetStructure().GetDimension()==3)
							PlotPhononFunction = & PlotPhononFrequencies_3DCubic;
						else if(solver[0]->GetStructure().GetDimension()==2)
							PlotPhononFunction = & PlotPhononFrequencies_2DCubic;
						else
							std::cerr<<"Unrecongnized symmetry!\n";
					}
					else if(strcmp(tempstring, "FCC")==0)
					{
						if(solver[0]->GetStructure().GetDimension()==3)
							PlotPhononFunction = & PlotPhononFrequencies_3DFCC;
						else
							std::cerr<<"Unrecongnized symmetry!\n";
					}
					else if(strcmp(tempstring, "Rectangle")==0 && solver[0]->GetStructure().GetDimension()==2)
					{
						PlotPhononFunction = & PlotPhononFrequencies_2DRectangle;
					}
					else if(strcmp(tempstring, "Tetragonal")==0 && solver[0]->GetStructure().GetDimension()==3)
					{
						PlotPhononFunction = & PlotPhononFrequencies_3DTetragonal;
					}
					else
						std::cerr<<"Unrecongnized symmetry!\n";
				}
			}

			else if(strcmp(tempstring, "Iterate")==0)
			{
				if(solver[0] == nullptr)
				{
					std::cerr<<"Specify Target before iteration!\n";
				}
				else if(pParam == nullptr)
				{
					std::cerr<<"Specify Parameterization before iteration!\n";
				}
				else if(NumInstances.size()==0)
				{
					std::cerr<<"Input some MC instances before iteration!\n";
				}
				else
				{
					if(PlotPhononFunction == nullptr)
					{
						std::cerr<<"Warning: No symmetry specified, will skip phonon frequency calculation!\n";
					}
					if(CoolingScheduleRescale==0.0)
						InvThermoIteration(solver, pParam, ReOptimizeTime, NumInstances, AnnealTargetNumRequirement, PlotPhononFunction, 1.0, &Pressure, RandomSeed, AllowSkipOptimization);
					else
						InvThermoIteration(solver, pParam, ReOptimizeTime, NumInstances, AnnealTargetNumRequirement, PlotPhononFunction, CoolingScheduleRescale, &Pressure, RandomSeed, AllowSkipOptimization);
				}
			}

			else if(strcmp(tempstring, "Verify")==0)
			{
				if(solver[0] == nullptr)
				{
					std::cerr<<"Specify Target before final check!\n";
					break;
				}
				if(pParam == nullptr)
				{
					std::cerr<<"Specify Parameterization before final check!\n";
					break;
				}
				if(PlotPhononFunction == nullptr)
				{
					std::cerr<<"Specify the Symmetry before final check!\n";
					break;
				}

				size_t NumP, NumC;
				std::cout<<"Input MC Num of Particles:";
				std::cin>>NumP;
				std::cout<<"Input MC Num of Configurations:";
				std::cin>>NumC;
				if(CoolingScheduleRescale==0.0)
					CheckResult(solver, pParam, NumP, NumC, PlotPhononFunction, 1.0);
				else
					CheckResult(solver, pParam, NumP, NumC, PlotPhononFunction, CoolingScheduleRescale);
			}
			else if(strcmp(tempstring, "SearchDegeneracy")==0)
			{
				if(solver[0] == nullptr)
				{
					std::cerr<<"Specify Target before final check!\n";
					break;
				}
				size_t TrialTime;
				std::cout<<"Input Trial Time:";
				std::cin>>TrialTime;
				size_t NumParticle;
				std::cout<<"Input Num of Particles:";
				std::cin>>NumParticle;
				std::vector<Configuration> results;
				RandomGenerator gen(RandomSeed);
				SearchDegenerate(solver[0], gen, TrialTime, results, NumParticle);
				Plot(std::string(solver[0]->Tag()), MultiplicateStructure(solver[0]->GetStructure(), 3));
				for(auto iter=results.begin(); iter!=results.end(); iter++)
				{
					std::string prefix(solver[0]->Tag());
					prefix+=std::string("_nbr_");
					prefix+=std::to_string(static_cast<long double>(NumParticle));
					prefix+=std::string("_degenerate_");
					prefix+=std::to_string(static_cast<long double>(iter-results.begin()));
					Plot(prefix, MultiplicateStructure(*iter, 3));
					Output(prefix, *iter);
				}
			}
			else if(strcmp(tempstring, "AllowSkipOptimization")==0)
			{
				std::cin>>AllowSkipOptimization;
			}
			else
			{
				std::cout<<"Unrecognized command!\n";
				std::cin.clear();
			}
			ifile.ignore(1000,'\n');
			if(ifile.eof()) return 2;
		}
	}
	catch (char const * a){std::cerr<<"error:"<<a<<std::endl;}
	return 1;
}


struct MC_Result
{
	double V, H;
	size_t Number;
};


const double m=0.0392304845;
const double sqrt3=std::sqrt(3.0);
class QuinticPotential : public IsotropicPairPotential
{
	double L, K, Frb;
	double F(double x)
	{
		double rmsqrt3=x-sqrt3;
		return (L*rmsqrt3*rmsqrt3+K*rmsqrt3-1)*rmsqrt3*rmsqrt3*rmsqrt3;
	}
public:
	QuinticPotential() : IsotropicPairPotential(2, sqrt3)
	{
		double sqrt3m1=sqrt3-1.0;
		this->L=4*m*std::pow(sqrt3m1, -5.0)-std::pow(sqrt3m1, -2.0);
		this->K=5*m*std::pow(sqrt3m1, -4.0)-2*std::pow(sqrt3m1, -1.0);
		double Rb=sqrt3-3.0*sqrt3m1/5.0/(1-4*m*std::pow(sqrt3m1, -3.0));
		this->Frb=this->F(Rb);
	}
	double IsotropicEnergyFormula(double distance)
	{
		return F(distance)/this->Frb;
	}
};
int OneDMC()
{
	const size_t NumParticle=1000;
	const size_t NumPlotParticle=50;
	const double Pressure=1.0;
	const double MinDistance=0.1;
	::AlsoWriteEPS=true;
	RandomGenerator gen;
	RnPotential pot(1, 1.0, -2, 100);
	GeometryVector basis=GeometryVector( (double)(NumParticle) );
	Configuration list(1, &basis, 1);
	for(size_t i=0; i<NumParticle; i++)
	{
		GeometryVector t(1);
		t.x[0]=((double)(i)+gen.RandomDouble())/NumParticle;
		list.Insert("A", t);
	}
	double E=pot.Energy(list);
	MonteCarlo<Configuration, Potential *> system(list, gen.RandomDouble()*65535, E);
	AtomMove ma;
	ma.MinDistance=MinDistance;

	for(double Temperature = 0.5; Temperature > 0.45; Temperature -=0.1)
	{
		std::cout<<"Temperature="<<Temperature<<'\n';
		for(size_t i=0; i<20; i++)
		{
			std::cout<<pot.Energy(list)<<'\n';
			system.Move(Temperature, 10000*NumParticle, ma, &pot);
		}

		std::vector<GeometryVector> results;
		std::string filename;

		IsotropicStructureFactor( [&list, &ma, &pot, &system, &NumParticle, &Temperature] (size_t i) ->Configuration
		{
			system.Move(Temperature, 100*NumParticle, ma, &pot);
			return list;
		}, 20, 2.0, 2.0, results);
		filename=std::to_string( (long double)(Temperature) )+std::string("_Sk.txt");
		std::fstream file1(filename.c_str(), std::fstream::out);
		for(auto iter=results.begin(); iter!=results.end(); iter++)
			file1<<(*iter);


		/*
		IsotropicTwoPairCorrelation( [&list, &ma, &pot, &system] (size_t i) ->Configuration
		{
		system.Move(Temperature, 100*NumParticle, ma, &pot);
		return list;
		}, 20, 5.0, results);
		filename=std::to_string( (long double)(Temperature) )+std::string("_g2.txt");
		std::fstream file2(filename.c_str(), std::fstream::out);
		for(auto iter=results.begin(); iter!=results.end(); iter++)
		file2<<(*iter);

		std::cout<<"outputting configuration\n";
		filename=std::to_string( (long double)(Temperature) )+std::string("_config");
		std::vector<GeometryVector> temp;
		temp.push_back(GeometryVector((double)(NumPlotParticle)/NumParticle, 0.0));
		temp.push_back(GeometryVector(0.0, (double)(NumPlotParticle)/NumParticle/20));
		Configuration templ(2, &temp[0], 1000);

		list.IterateThroughNeighbors(GeometryVector((double)(NumPlotParticle)/NumParticle/2.0), NumPlotParticle/2.0, [&templ](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
		{
		GeometryVector t(2);
		t.x[0]=shift.x[0]/NumPlotParticle;
		t.x[1]=0.5;
		templ.Insert("A", t);
		} );
		::Plot(filename, templ);
		*/
	}
	return 0;
}
double LowestEigenValueOfDeformations(Configuration list, Potential & pot, double Pressure);
//int test2()
//{
//	std::fstream paramfile("CheckResult_Parameters.txt", std::fstream::in);
//	LoadStructure sol("2D/AntiKagome.pos");
//	ParameterSet * p1=new withGaussianInvParameters(&sol);
//	p1->Evolve();
//	ParameterSet * p2=new SoftCoreParameterization(&sol);
//	CombinedParameterization param(p1, p2);
//	param.Read(paramfile);
//	IsotropicPairPotential * pPot=param.GetPotential();
//	pPot->IsotropicPairEnergy(2.0);
//	delete pPot;
//	return 0;
//}
int PlotAll()
{
	for(size_t i =0; ; i++)
	{
		std::stringstream stemp;
		char PrefixBuffer[64];
		stemp<<"CheckResult_MCAnneal_"<<i<<".pos";
		stemp.getline(PrefixBuffer, 63);
		std::fstream infile(PrefixBuffer, std::fstream::in);
		if(infile.good()==false)
			break;
		Configuration a=ReadPos(infile);
		::Plot(std::string(PrefixBuffer), a);
		infile.close();
	}
	return 0;
}
int RescaleAnnealings()
{
	const double factor=0.5;
	for(size_t i =0; ; i++)
	{
		std::stringstream stemp;
		char PrefixBuffer[64];
		stemp<<"Annealing_"<<i<<".pos";
		stemp.getline(PrefixBuffer, 63);
		std::fstream infile(PrefixBuffer, std::fstream::in);
		if(infile.good()==false)
			break;
		Configuration a=ReadPos(infile);
		{
			std::stringstream stemp;
			stemp<<"BackUp_Annealing_"<<i;
			stemp.getline(PrefixBuffer, 63);
			::Output(PrefixBuffer, a);
		}
		std::vector<GeometryVector> bas;
		for(DimensionType j=0; j<a.GetDimension(); j++)
			bas.push_back(a.GetBasisVector(j)*factor);
		a.ChangeBasisVector(&bas[0]);
		{
			std::stringstream stemp;
			stemp<<"Annealing_"<<i;
			stemp.getline(PrefixBuffer, 63);
			::Output(PrefixBuffer, a);
		}
		infile.close();
	}
	return 0;
}
int test()
{
	std::vector<GeometryVector> bas;
	bas.push_back(GeometryVector(1.0, 0.0));
	bas.push_back(GeometryVector(0.0, 1.0));
	PeriodicCellList<Empty> lattice(2, &bas[0], 1.0, false);
	lattice.Insert(Empty(), GeometryVector(0.0, 0.0));

	const double Radius1=30;
	const double Radius2=2;
	for(int i=0; i<10; i++)
	{
		GeometryVector center(Radius1*std::cos(2*pi*i/10.0), Radius1*std::sin(2*pi*i/10.0));
		lattice.IterateThroughNeighbors(center, Radius2, [&center](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t source) ->void
		{
			GeometryVector t=center+shift;
			std::cout<<t.x[0]<<" \t"<<t.x[1]<<'\n';
		});
	}
	return 0;
}
int test3()
{
	withGaussianInvPotential pot(2, 3.0);
	pot.invterms.push_back(InvTerm(1.0, 12.0));
	pot.invterms.push_back(InvTerm(-2.0, 6.0));
	pot.wells.push_back(GaussianWell(0.7, 0.2, std::sqrt(3.0)));

	std::vector<LatticeSumSolver *> solver;
	solver.push_back(new LoadStructure("2D/Square.pos"));
	solver.push_back(new LoadStructure("2D/Rectangle1.1.pos"));
	solver.push_back(new LoadStructure("2D/Triangle.pos"));
	solver.push_back(new LoadStructure("2D/Graphene.pos"));
	solver.push_back(new LoadStructure("2D/Kagome.pos"));

	::AlsoWriteEPS=true;
	::Plots(&pot, &solver, "test", 1);
	return 0;
}
int test4()
{
	std::vector<LatticeSumSolver *> solver;
	ParameterSet * pParam = nullptr;

	solver.push_back(new LoadStructure("2D/Kagome.pos"));
	solver.push_back(new LoadStructure("2D/Triangle.pos"));

	pParam = new GaussianCoreParameterization(solver[0]);

	size_t NumParam, NumPosOutput;
	RestoreStatus(*pParam, solver, NumParam, NumPosOutput);
	IsotropicPairPotential * pot = pParam->GetPotential();
	Configuration tar=solver[0]->GetStructure();
	double Pressure=GetPressure(solver[0], tar.PeriodicVolume()/tar.NumParticle(), *pot);
	RelaxStructure(tar, *pot, Pressure, 0.0);
	PlotPhononFrequencies_2DHexagonal(tar, *pot, "test", 1.0, 100);
	std::cout<<"Deigen:"<<LowestEigenValueOfDeformations(tar, *pot, Pressure)<<'\n';
	std::cout<<"Ceigen:"<<LowestEigenValueOfStiffness(tar, *pot, Pressure)<<'\n';
	std::cout<<"Ediff:"<<EnthalpyDifference::EnthalpyDifference(*pParam, solver, tar.PeriodicVolume()/tar.NumParticle(), Pressure);
	return 0;
}
int test5()
{
	size_t NumParticle1=3;
	size_t NumParticle2=3;
	RandomGenerator gen;
	std::vector<std::pair<Configuration, Configuration> > result;
	SearchDegenerate2(gen, 20, result, NumParticle1, NumParticle2, 2);

	double tempRadius=20;
	std::swap(tempRadius, ::TwoBodyDistance_MaxLength);

	for(auto iter=result.begin(); iter!=result.end(); iter++)
	{
		std::string prefix("");
		prefix+=std::string("nbr_");
		prefix+=std::to_string(static_cast<long double>(NumParticle1));
		prefix+=std::string("_nbr_");
		prefix+=std::to_string(static_cast<long double>(NumParticle2));
		prefix+=std::string("_degenerate_");
		prefix+=std::to_string(static_cast<long double>(iter-result.begin()));
		FromStructure sa(iter->first);
		FromStructure sb(iter->second);

		SolverDistance(&sa, &sb);

		std::string namea=prefix+std::string("_a");
		Plot(namea, MultiplicateStructure(iter->first, 3));
		Output(namea, iter->first);
		std::fstream filea((namea+std::string(".ThetaSeries")).c_str(), std::fstream::out);
		OutputThetaSeries(sa, filea);

		std::string nameb=prefix+std::string("_b");
		Plot(nameb, MultiplicateStructure(iter->second, 3));
		Output(nameb, iter->second);
		std::fstream fileb((nameb+std::string(".ThetaSeries")).c_str(), std::fstream::out);
		OutputThetaSeries(sb, fileb);

	}
	std::swap(tempRadius, ::TwoBodyDistance_MaxLength);
	return 0;
}


////potential corresponding to Optimization_37 when optimizing for H
//class CaF2Potential : public IsotropicPairPotential
//{
//public:
//	CaF2Potential() : IsotropicPairPotential(3, 2.0564)
//	{
//	}
//	virtual double IsotropicEnergyFormula(double x)
//	{        
//		double V=2.9340e-03*std::pow((1/x), 12)+8.3963e-01+ 3.6976e-01*x-1.3150e-01*x*x- 2.1869e-03*x*x*x+1.5010e-03*std::pow(x, 4);
//		V=V*(x-2.0564)*(x-2.0564);
//		V=V*exp( -1.8682e-01*x*x );
//		return V;
//	}
//};
//class Rectangle_3_Potential : public IsotropicPairPotential
//{
//public:
//	Rectangle_3_Potential() : IsotropicPairPotential(2, 4)
//	{
//	}
//	virtual double IsotropicEnergyFormula(double x)
//	{        
//		double V=2.1639e-02*std::pow((1/x), 12)- 2.6107e-01+ 3.1488e-01*x;
//		V=V*(x-6.4)*(x-6.4);
//		V=V*exp( -0.788571*x*x );
//		return V;
//	}
//};
//class Kagome_Potential : public IsotropicPairPotential
//{
//public:
//	Kagome_Potential() : IsotropicPairPotential(2, 2.03643)
//	{
//	}
//	virtual double IsotropicEnergyFormula(double x)
//	{        
//		double V=5.9860e-02*std::pow((1/x), 12)-1.2811+ 2.1521*x;
//		V=V*(x-2.03643)*(x-2.03643);
//		return V;
//	}
//};
//
//class Rectangle_7_Potential : public IsotropicPairPotential
//{
//public:
//	Rectangle_7_Potential() : IsotropicPairPotential(2, 2.25245)
//	{
//	}
//	virtual double IsotropicEnergyFormula(double x)
//	{        
//		double V=3.0058e-03*std::pow((1/x), 12)+ 6.9293e-01-3.0361e-01*x+9.3960e-02*x*x-3.6154e-01*x*x*x+ 8.2231e-01*x*x*x*x+4.3741e-02*x*x*x*x*x;
//		V=V*(x-2.25245)*(x-2.25245);
//		V=V*exp( -0.440952*x*x );
//		return V;
//	}
//};
//class Rectanglepi_Potential : public IsotropicPairPotential
//{
//public:
//	Rectanglepi_Potential() : IsotropicPairPotential(2, 3.41025)
//	{
//	}
//	virtual double IsotropicEnergyFormula(double x)
//	{        
//		double V= 1.1416e-02*std::pow((1/x), 12)-1.1117+3.3164*x-3.1330*x*x+1.2578*x*x*x-1.6340e-01*x*x*x*x;
//		V=V*(x-3.41025)*(x-3.41025);
//		V=V*exp( -0.0309012 *x*x );
//		return V;
//	}
//};
//class RectangularKagome_Potential : public IsotropicPairPotential
//{
//public:
//	RectangularKagome_Potential() : IsotropicPairPotential(2, 3.050295)
//	{
//	}
//	virtual double IsotropicEnergyFormula(double x)
//	{        
//		double V=(0.012352*x+0.27370)*(x-3.050295)*(x-3.050295)*exp( -0.086364*x*x );
//
//		double temp=(x-0.99953)/0.024893;
//		V-= 0.092965*std::exp((-1.0)*temp*temp);
//		V+= 3.8032e-04*std::pow((1/x), 12) - 1.0430e-02*std::pow((1/x), 6);
//		V+=1.2956e-05;
//
//		return V;
//	}
//};
////potential corresponding to Optimization_37 when optimizing for H
//class CaF2_P_Potential : public IsotropicPairPotential
//{
//public:
//	CaF2_P_Potential() : IsotropicPairPotential(3, 0.890432)
//	{
//	}
//	virtual double IsotropicEnergyFormula(double x)
//	{        
//		double V=std::pow((0.34641/x), 12)+19.6654+20*x-16.4269*x*x-0.630861*x*x*x+std::pow(x, 4.0);
//		V=V*(x-0.890432)*(x-0.890432);
//		V=V*std::exp( -0.996388*x*x );
//		return V;
//	}
//};
//
//#include <nlopt.h>
//namespace AspectRatioAtDifferentPressure
//{
//	struct optimize_temp
//	{
//		IsotropicPairPotential * pPot;
//		double * pPressure;
//	};
//	double Obj(unsigned int n, const double *x, double *grad, void *aux)
//	{
//		assert(n==2);
//		std::vector<GeometryVector> bases;
//		bases.push_back(GeometryVector(x[0], 0.0));
//		bases.push_back(GeometryVector(0.0, x[1]));
//		Configuration rect(2, &bases[0],  1.0, false);
//		rect.Insert("A", GeometryVector(0.0, 0.0));
//		optimize_temp * ptemp = reinterpret_cast<optimize_temp *>(aux);
//		return ptemp->pPot->Energy(rect)+(*ptemp->pPressure)*x[0]*x[1];
//	}
//	Configuration GetOptimizedRectangle(IsotropicPairPotential * pPot, double Pressure, double & AspectRatio)
//	{
//		optimize_temp aux;
//		aux.pPot=pPot;
//		aux.pPressure=&Pressure;
//		nlopt_opt opt;
//		opt=nlopt_create(NLOPT_LN_NELDERMEAD, 2);
//		nlopt_set_ftol_abs(opt, 1e-8);
//		nlopt_set_min_objective(opt, Obj, reinterpret_cast<void *>(&aux));
//		double sidelengths[2] = {2.0, 1.0};
//		double lbound[2]={1.0, 0.5};
//		double ubound[2]={5.0, 2.0};
//		nlopt_set_upper_bounds(opt, ubound);
//		nlopt_set_lower_bounds(opt, lbound);
//		double minf;
//		nlopt_optimize(opt, sidelengths, &minf);
//		AspectRatio=sidelengths[0]/sidelengths[1];
//		std::vector<GeometryVector> bases;
//		bases.push_back(GeometryVector(sidelengths[0], 0.0));
//		bases.push_back(GeometryVector(0.0, sidelengths[1]));
//		Configuration rect(2, &bases[0],  1.0, false);
//		rect.Insert("A", GeometryVector(0.0, 0.0));
//		return rect;
//	}
//	int test6()
//	{
//		//Rectangle2_Potential Pot;
//		const size_t NumConfig=3;
//		const DimensionType dim=2;
//		const size_t MCNumParticle=36;
//		//const double MCPressure=335.107;//CaF2_H
//		//const double MCPressure=25.812025934808;//Rectangle_7
//		//const double MCPressure=6.2599488453739;//Rectanglepi
//		//const double MCPressure=5.75544;//Rectangle_3
//		const double CoolSchedule=1.5;
//		const double MinDistance=0.5;
//		AlsoWriteEPS=true;
//
//		LoadStructure target("2D/Rectangle2.pos");
//		//IsotropicPairPotential * pPot0 = new Rectangle_7_Potential();
//		//CorrectedIsotropicPairPotential pot(pPot0);
//		Rectangle_7_Potential pot;
//		IsotropicPairPotential * pPot=&pot;
//		IsotropicPairPotential & Pot = (*pPot);
//
//		Configuration tutu=target.GetStructure();
//		double MCPressure = GetPressure(&target, tutu.PeriodicVolume()/tutu.NumParticle(), Pot);
//		//std::cout<<(Pot.Energy(tutu)+MCPressure*tutu.PeriodicVolume())/tutu.NumParticle()<<'\n';
//
//		std::vector<double> pressures;
//		for(double t=0.82; t<2.22; t*=1.005)
//			pressures.push_back(t);
//		std::vector<Configuration> results;
//		results.resize(pressures.size(), tutu);
//		/*
//		LoadStructure tri("2D/Triangle.pos");
//		for(double i=0.3; i<5; i*=1.01)
//		std::cout<<i<<" \t"<<target.LatticeSum(i, Pot)+i*MCPressure<<'\n';
//		for(double i=0.3; i<3.1; i*=1.001)
//		std::cout<<i<<" \t"<<Pot.SecondDerivative(GeometryVector(i, 0.0), 0, 0)<<'\n';
//		std::cout<<std::sqrt(5.0)<<" \t"<<Pot.SecondDerivative(GeometryVector(std::sqrt(5.0), 0.0), 0, 0)<<'\n';
//		return 0;
//		*/
//
//		signed int end=pressures.size();
//		std::cout<<"Num Configurations:"<<end<<'\n';
//#pragma omp parallel for num_threads(7)
//		for(signed int i=0; i<end; i++)
//		{
//			ParticleSimulatedAnnealing a(dim, MCNumParticle, 1.0, i, pressures[i], Pot, i==0?&std::cout:nullptr, CoolSchedule, MinDistance);
//			a.Anneal(Pot);
//			std::stringstream Number;
//			Number<<i;
//			std::string name("CheckResult_MCAnneal_");
//			std::string name2("");
//			Number>>name2;
//			name+=name2;
//			Output(name.c_str(), *a.pList);
//			std::cout<<i<<"th result completed! E/N="<<Pot.Energy(*a.pList)/MCNumParticle<<", V/N="<<a.pList->PeriodicVolume()/MCNumParticle<<", H/N="<<(Pot.Energy(*a.pList)+MCPressure*a.pList->PeriodicVolume())/MCNumParticle<<'\n';
//			logfile<<i<<"th result completed! E/N="<<Pot.Energy(*a.pList)/MCNumParticle<<", V/N="<<a.pList->PeriodicVolume()/MCNumParticle<<", H/N="<<(Pot.Energy(*a.pList)+MCPressure*a.pList->PeriodicVolume())/MCNumParticle<<'\n';
//
//			std::cout<<"Calculate Theta Series\n";
//			logfile<<"Calculate Theta Series\n";
//			LoadStructure l((name+std::string(".pos")).c_str());
//			LJPotential lj(dim, 1, 1, 4);
//			l.LatticeSum(1.0, lj);
//			std::fstream ofile((name+std::string(".ThetaSeries")).c_str(), std::fstream::out);
//			OutputThetaSeries(l, ofile);
//
//			Plot(name, *a.pList);
//			results[i]=*a.pList;
//		}
//		for(signed int i=0; i<end; i++)
//		{
//			FromStructure temp(results[i]);
//			double AspectRatio;
//			Configuration rect=GetOptimizedRectangle(pPot, pressures[i], AspectRatio);
//			FromStructure std(rect);
//			std::cout<<pressures[i]<<" \t"<<rect.GetBasisVector(0).x[0]<<" \t"<<rect.GetBasisVector(1).x[1]<<" \t"<<SolverDistance(&std, &temp)<<'\n';
//		}
//
//		return 0;
//	}
//};
//namespace PoissonRatio
//{
//	struct optimize_temp
//	{
//		IsotropicPairPotential * pPot;
//		double * pPressure;
//	};
//	double Obj(unsigned int n, const double *x, double *grad, void *aux)
//	{
//		assert(n==2);
//		std::vector<GeometryVector> bases;
//		bases.push_back(GeometryVector(x[0], 0.0));
//		bases.push_back(GeometryVector(0.0, x[1]));
//		Configuration rect(2, &bases[0],  1.0, false);
//		rect.Insert("A", GeometryVector(0.0, 0.0));
//		optimize_temp * ptemp = reinterpret_cast<optimize_temp *>(aux);
//		return ptemp->pPot->Energy(rect)+(*ptemp->pPressure)*x[0]*x[1];
//	}
//	double GetOptimizedRectangle(IsotropicPairPotential * pPot, double Pressure, double & a)
//	{
//		optimize_temp aux;
//		aux.pPot=pPot;
//		aux.pPressure=&Pressure;
//		nlopt_opt opt;
//		opt=nlopt_create(NLOPT_LN_NELDERMEAD, 2);
//		nlopt_set_ftol_abs(opt, 1e-11);
//		nlopt_set_min_objective(opt, Obj, reinterpret_cast<void *>(&aux));
//		double sidelengths[2] = {2.0, 1.0};
//		double lbound[2]={1.0, 0.5};
//		double ubound[2]={5.0, 2.0};
//		sidelengths[0]=a;
//		lbound[0]=a;
//		ubound[0]=a;
//		nlopt_set_upper_bounds(opt, ubound);
//		nlopt_set_lower_bounds(opt, lbound);
//		double minf;
//		nlopt_optimize(opt, sidelengths, &minf);
//		return sidelengths[1];
//	}
//	int test6()
//	{
//		const double pressure=1.12901;
//		LoadStructure target("2D/Rectangle2.pos");
//		//IsotropicPairPotential * pPot0 = new Rectangle_7_Potential();
//		//CorrectedIsotropicPairPotential pot(pPot0);
//		Rectangle_7_Potential pot;
//		IsotropicPairPotential * pPot=&pot;
//		IsotropicPairPotential & Pot = (*pPot);
//
//		Configuration tutu=target.GetStructure();
//		double MCPressure = GetPressure(&target, tutu.PeriodicVolume()/tutu.NumParticle(), Pot);
//		//std::cout<<(Pot.Energy(tutu)+MCPressure*tutu.PeriodicVolume())/tutu.NumParticle()<<'\n';
//
//		for(double a=1.999; a<2.001; a+=0.0001)
//		{
//			double b=GetOptimizedRectangle(pPot, pressure, a);
//			std::cout<<a<<" \t"<<b<<'\n';
//		}
//
//		return 0;
//	}
//};
//int test7()
//{
//	double LowestEigenValueOfDeformations(Configuration list, Potential & pot, double Pressure);
//	void DeformationsMatrixPrincipleMinorDeterminants(Configuration list, Potential & pot, double Pressure, double * results);
//	LoadStructure target("2D/Kagome.pos");
//	GaussianCoreParameterization para(&target);
//	std::fstream parafile("Optimization_17_Parameters.txt", std::fstream::in);
//	para.Read(parafile);
//	IsotropicPairPotential * pPot=para.GetPotential();
//	Configuration con=target.GetStructure();
//	double Pressure=GetPressure(&target, con.PeriodicVolume()/con.NumParticle(), *pPot);
//	double res[100];
//	LowestEigenValueOfDeformations(con, *pPot, Pressure);
//	DeformationsMatrixPrincipleMinorDeterminants(con, *pPot, Pressure, res);
//	return 0;
//}
//
//int test8()
//{
//	//::AlsoWriteEPS=true;
//	{
//		CaF2Potential pot;
//		std::vector<LatticeSumSolver *> solvers;
//		solvers.push_back(new LoadStructure("3D/CaF2_rescaled.pos"));
//		::Plots(&pot, &solvers, "CaF2_H_");
//		PlotPhononFrequencies_3DFCC(solvers[0]->GetStructure(), pot, "CaF2_H_", 1.0, 500);
//		Configuration tar=solvers[0]->GetStructure();
//		double vpn=tar.PeriodicVolume()/tar.NumParticle();
//		double Pressure=GetPressure(solvers[0], vpn, pot);
//		std::fstream outfile("CaF2_ElasticConstants.txt", std::fstream::out);
//		PrintAllElasticConstants(outfile, solvers[0]->GetStructure(), pot, Pressure);
//		outfile<<"Pressure="<<Pressure<<'\n';
//	}
//	{
//		double temp=200;
//		std::swap(temp, ::PlotPhononFrequencies_MaxFrequencySquared_Override);
//		Rectangle_3_Potential pot;
//		std::vector<LatticeSumSolver *> solvers;
//		solvers.push_back(new LoadStructure("2D/Rectangle2.pos"));
//		::Plots(&pot, &solvers, "Rectangle_H_3");
//		PlotPhononFrequencies_2DRectangle(solvers[0]->GetStructure(), pot, "Rectangle_H_3", 1.0, 500);
//		std::swap(temp, ::PlotPhononFrequencies_MaxFrequencySquared_Override);
//		Configuration tar=solvers[0]->GetStructure();
//		double vpn=tar.PeriodicVolume()/tar.NumParticle();
//		double Pressure=GetPressure(solvers[0], vpn, pot);
//		std::fstream outfile("Rectangle_H_3_ElasticConstants.txt", std::fstream::out);
//		PrintAllElasticConstants(outfile, solvers[0]->GetStructure(), pot, Pressure);
//		outfile<<"Pressure="<<Pressure<<'\n';
//	}
//	{
//		Rectangle_7_Potential pot;
//		std::vector<LatticeSumSolver *> solvers;
//		solvers.push_back(new LoadStructure("2D/Rectangle2.pos"));
//		::Plots(&pot, &solvers, "Rectangle_H_7");
//		PlotPhononFrequencies_2DRectangle(solvers[0]->GetStructure(), pot, "Rectangle_H_7", 1.0, 500);
//		Configuration tar=solvers[0]->GetStructure();
//		double vpn=tar.PeriodicVolume()/tar.NumParticle();
//		double Pressure=GetPressure(solvers[0], vpn, pot);
//		std::fstream outfile("Rectangle_H_7_ElasticConstants.txt", std::fstream::out);
//		PrintAllElasticConstants(outfile, solvers[0]->GetStructure(), pot, Pressure);
//		outfile<<"Pressure="<<Pressure<<'\n';
//	}
//	{
//		Rectanglepi_Potential pot;
//		std::vector<LatticeSumSolver *> solvers;
//		solvers.push_back(new LoadStructure("2D/Rectanglepi.pos"));
//		::Plots(&pot, &solvers, "Rectanglepi_H");
//		PlotPhononFrequencies_2DRectangle(solvers[0]->GetStructure(), pot, "Rectanglepi_H", 1.0, 500);
//		Configuration tar=solvers[0]->GetStructure();
//		double vpn=tar.PeriodicVolume()/tar.NumParticle();
//		double Pressure=GetPressure(solvers[0], vpn, pot);
//		std::fstream outfile("Rectanglepi_H_ElasticConstants.txt", std::fstream::out);
//		PrintAllElasticConstants(outfile, solvers[0]->GetStructure(), pot, Pressure);
//		outfile<<"Pressure="<<Pressure<<'\n';
//	}
//	{
//		Kagome_Potential pot;
//		std::vector<LatticeSumSolver *> solvers;
//		solvers.push_back(new LoadStructure("2D/Kagome.pos"));
//		::Plots(&pot, &solvers, "Kagome_H");
//		PlotPhononFrequencies_2DHexagonal(solvers[0]->GetStructure(), pot, "Kagome_H", 1.0, 500);
//		Configuration tar=solvers[0]->GetStructure();
//		double vpn=tar.PeriodicVolume()/tar.NumParticle();
//		double Pressure=GetPressure(solvers[0], vpn, pot);
//		std::fstream outfile("Kagome_ElasticConstants.txt", std::fstream::out);
//		PrintAllElasticConstants(outfile, solvers[0]->GetStructure(), pot, Pressure);
//		outfile<<"Pressure="<<Pressure<<'\n';
//	}
//
//	{
//		double temp=2000;
//		std::swap(temp, ::PlotPhononFrequencies_MaxFrequencySquared_Override);
//		RectangularKagome_Potential pot;
//		std::vector<LatticeSumSolver *> solvers;
//		solvers.push_back(new LoadStructure("2D/RectangularKagome.pos"));
//		::Plots(&pot, &solvers, "RectangularKagome_H");
//		PlotPhononFrequencies_2DRectangle(solvers[0]->GetStructure(), pot, "RectangularKagome_H", 1.0, 500);
//		std::swap(temp, ::PlotPhononFrequencies_MaxFrequencySquared_Override);
//		Configuration tar=solvers[0]->GetStructure();
//		double vpn=tar.PeriodicVolume()/tar.NumParticle();
//		double Pressure=GetPressure(solvers[0], vpn, pot);
//		std::fstream outfile("RectangularKagome_H_ElasticConstants.txt", std::fstream::out);
//		PrintAllElasticConstants(outfile, solvers[0]->GetStructure(), pot, Pressure);
//		outfile<<"Pressure="<<Pressure<<'\n';
//	}
//
//
//	return 0;
//}
//int GenerateCaF2InitConfiguration()
//{
//	//std::ostream & out=std::cout;
//	std::fstream out("coord.txt", std::fstream::out);
//	out.precision(14);
//	Configuration caf2=ReadPos("3D/CaF2");
//	size_t count=0;
//	GeometryVector basis[3]={GeometryVector(10.0, 0.0, 0.0), GeometryVector(0.0, 10.0, 0.0), GeometryVector(0.0, 0.0, 10.0)};
//	Configuration init(3, basis, 1.0, true);
//	caf2.IterateThroughNeighbors(GeometryVector(3), std::sqrt(51.1), [&out, &init, &count](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void
//	{
//		if(shift.x[0]>(-0.0001) && shift.x[0]<0.9999)
//			if(shift.x[1]>(-5.0001) && shift.x[1]<4.9999)
//				if(shift.x[2]>(-5.0001) && shift.x[2]<4.9999)
//				{
//					init.Insert("B", shift*0.1);
//					out<<shift.x[0]<<" \t"<<shift.x[1]<<" \t"<<shift.x[2]<<'\n';
//					count++;
//				}
//	});
//	assert(count==1200);
//	RandomGenerator gen(0);
//	for(; count<12000; )
//	{
//		double x=(gen.RandomDouble()-0.5)*10.0;
//		if(x>(-0.0001) && x<0.99999)
//			continue;
//		GeometryVector t(x*0.1, gen.RandomDouble()-0.5, gen.RandomDouble()-0.5);
//		bool HasNeighbor=false;
//		init.IterateThroughNeighbors(t, 0.3, [&HasNeighbor](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void
//		{
//			HasNeighbor=true;
//		});
//		if(HasNeighbor==true)
//			continue;
//
//		init.Insert("A", t);
//		GeometryVector t2=t*10;
//		out<<t2.x[0]<<" \t"<<t2.x[1]<<" \t"<<t2.x[2]<<" \n";
//		count++;
//	}
//	return 0;
//}
//
//#include "DisplaySpheres.h"
//line GetLine(double x1, double y1, double z1, double x2, double y2, double z2)
//{
//	line result;
//	result.x1=x1;
//	result.y1=y1;
//	result.z1=z1;
//	result.x2=x2;
//	result.y2=y2;
//	result.z2=z2;
//	result.blue=0.0;
//	result.red=0.0;
//	result.green=0.0;
//	result.transparency=1.0;
//	result.width=0.05;
//	return result;
//}
//sphere GetCa(double x, double y, double z)
//{
//	sphere result;
//	result.radius=0.07*1.26;
//	result.x=x;
//	result.y=y;
//	result.z=z;
//	result.transparency=1.0;
//	result.blue=0.8;
//	result.red=0.0;
//	result.green=0.0;
//	return result;
//}
//sphere GetF(double x, double y, double z)
//{
//	sphere result;
//	result.radius=0.07*1.17;
//	result.x=x;
//	result.y=y;
//	result.z=z;
//	result.transparency=1.0;
//	result.blue=0.3;
//	result.red=1.0;
//	result.green=1.0;
//	return result;
//}
//sphere Geta(double x, double y, double z)
//{
//	sphere result;
//	result.radius=0.07;
//	result.x=x;
//	result.y=y;
//	result.z=z;
//	result.transparency=1.0;
//	result.blue=0.1;
//	result.red=1.0;
//	result.green=0.1;
//	return result;
//}
//sphere GetGray(double x, double y, double z)
//{
//	sphere result;
//	result.radius=0.07;
//	result.x=x;
//	result.y=y;
//	result.z=z;
//	result.transparency=1.0;
//	result.blue=0.1;
//	result.red=0.1;
//	result.green=0.1;
//	return result;
//}
//int DisplayCaF2()
//{
//	std::vector<line> Lines;
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 1.0, 0.0, 0.0));
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 0.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 0.0, 0.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 0.0, 1.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 1.0, 0.0, 1.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 1.0, 0.0, 0.0, 1.0, 1.0));
//	Lines.push_back(GetLine(0.0, 0.0, 1.0, 0.0, 1.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 0.0, 1.0, 0.0, 1.0));
//	Lines.push_back(GetLine(0.0, 0.0, 1.0, 1.0, 0.0, 1.0));
//	Lines.push_back(GetLine(1.0, 1.0, 0.0, 1.0, 1.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 1.0, 1.0, 1.0, 1.0));
//	Lines.push_back(GetLine(0.0, 1.0, 1.0, 1.0, 1.0, 1.0));
//	std::vector<sphere> spheres;
//	spheres.push_back(GetCa(0.0, 0.0, 0.0));
//	spheres.push_back(GetCa(1.0, 0.0, 0.0));
//	spheres.push_back(GetCa(0.0, 1.0, 0.0));
//	spheres.push_back(GetCa(0.0, 0.0, 1.0));
//	spheres.push_back(GetCa(1.0, 1.0, 0.0));
//	spheres.push_back(GetCa(1.0, 0.0, 1.0));
//	spheres.push_back(GetCa(0.0, 1.0, 1.0));
//	spheres.push_back(GetCa(1.0, 1.0, 1.0));
//	spheres.push_back(GetCa(0.5, 0.5, 0.0));
//	spheres.push_back(GetCa(0.5, 0.0, 0.5));
//	spheres.push_back(GetCa(0.0, 0.5, 0.5));
//	spheres.push_back(GetCa(0.5, 0.5, 1.0));
//	spheres.push_back(GetCa(0.5, 1.0, 0.5));
//	spheres.push_back(GetCa(1.0, 0.5, 0.5));
//	spheres.push_back(GetF(0.25, 0.25, 0.25));
//	spheres.push_back(GetF(0.75, 0.25, 0.25));
//	spheres.push_back(GetF(0.25, 0.75, 0.25));
//	spheres.push_back(GetF(0.25, 0.25, 0.75));
//	spheres.push_back(GetF(0.75, 0.75, 0.25));
//	spheres.push_back(GetF(0.75, 0.25, 0.75));
//	spheres.push_back(GetF(0.25, 0.75, 0.75));
//	spheres.push_back(GetF(0.75, 0.75, 0.75));
//	DisplaySpheres(spheres, Lines);
//	return 0;
//}
//int DisplayCaF2Single()
//{
//	std::vector<line> Lines;
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 1.0, 0.0, 0.0));
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 0.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 0.0, 0.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 0.0, 1.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 1.0, 0.0, 1.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 1.0, 0.0, 0.0, 1.0, 1.0));
//	Lines.push_back(GetLine(0.0, 0.0, 1.0, 0.0, 1.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 0.0, 1.0, 0.0, 1.0));
//	Lines.push_back(GetLine(0.0, 0.0, 1.0, 1.0, 0.0, 1.0));
//	Lines.push_back(GetLine(1.0, 1.0, 0.0, 1.0, 1.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 1.0, 1.0, 1.0, 1.0));
//	Lines.push_back(GetLine(0.0, 1.0, 1.0, 1.0, 1.0, 1.0));
//	std::vector<sphere> spheres;
//	spheres.push_back(Geta(0.0, 0.0, 0.0));
//	spheres.push_back(Geta(1.0, 0.0, 0.0));
//	spheres.push_back(Geta(0.0, 1.0, 0.0));
//	spheres.push_back(Geta(0.0, 0.0, 1.0));
//	spheres.push_back(Geta(1.0, 1.0, 0.0));
//	spheres.push_back(Geta(1.0, 0.0, 1.0));
//	spheres.push_back(Geta(0.0, 1.0, 1.0));
//	spheres.push_back(Geta(1.0, 1.0, 1.0));
//	spheres.push_back(Geta(0.5, 0.5, 0.0));
//	spheres.push_back(Geta(0.5, 0.0, 0.5));
//	spheres.push_back(Geta(0.0, 0.5, 0.5));
//	spheres.push_back(Geta(0.5, 0.5, 1.0));
//	spheres.push_back(Geta(0.5, 1.0, 0.5));
//	spheres.push_back(Geta(1.0, 0.5, 0.5));
//	spheres.push_back(Geta(0.25, 0.25, 0.25));
//	spheres.push_back(Geta(0.75, 0.25, 0.25));
//	spheres.push_back(Geta(0.25, 0.75, 0.25));
//	spheres.push_back(Geta(0.25, 0.25, 0.75));
//	spheres.push_back(Geta(0.75, 0.75, 0.25));
//	spheres.push_back(Geta(0.75, 0.25, 0.75));
//	spheres.push_back(Geta(0.25, 0.75, 0.75));
//	spheres.push_back(Geta(0.75, 0.75, 0.75));
//	DisplaySpheres(spheres, Lines);
//	return 0;
//}
//int DisplayDiamond()
//{
//	std::vector<line> Lines;
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 1.0, 0.0, 0.0));
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 0.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 0.0, 0.0, 0.0, 0.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 0.0, 1.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 1.0, 0.0, 1.0, 1.0, 0.0));
//	Lines.push_back(GetLine(0.0, 1.0, 0.0, 0.0, 1.0, 1.0));
//	Lines.push_back(GetLine(0.0, 0.0, 1.0, 0.0, 1.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 0.0, 1.0, 0.0, 1.0));
//	Lines.push_back(GetLine(0.0, 0.0, 1.0, 1.0, 0.0, 1.0));
//	Lines.push_back(GetLine(1.0, 1.0, 0.0, 1.0, 1.0, 1.0));
//	Lines.push_back(GetLine(1.0, 0.0, 1.0, 1.0, 1.0, 1.0));
//	Lines.push_back(GetLine(0.0, 1.0, 1.0, 1.0, 1.0, 1.0));
//	std::vector<sphere> spheres;
//	spheres.push_back(GetGray(0.0, 0.0, 0.0));
//	spheres.push_back(GetGray(1.0, 0.0, 0.0));
//	spheres.push_back(GetGray(0.0, 1.0, 0.0));
//	spheres.push_back(GetGray(0.0, 0.0, 1.0));
//	spheres.push_back(GetGray(1.0, 1.0, 0.0));
//	spheres.push_back(GetGray(1.0, 0.0, 1.0));
//	spheres.push_back(GetGray(0.0, 1.0, 1.0));
//	spheres.push_back(GetGray(1.0, 1.0, 1.0));
//	spheres.push_back(GetGray(0.5, 0.5, 0.0));
//	spheres.push_back(GetGray(0.5, 0.0, 0.5));
//	spheres.push_back(GetGray(0.0, 0.5, 0.5));
//	spheres.push_back(GetGray(0.5, 0.5, 1.0));
//	spheres.push_back(GetGray(0.5, 1.0, 0.5));
//	spheres.push_back(GetGray(1.0, 0.5, 0.5));
//	spheres.push_back(GetGray(0.25, 0.25, 0.25));
//	spheres.push_back(GetGray(0.75, 0.75, 0.25));
//	spheres.push_back(GetGray(0.75, 0.25, 0.75));
//	spheres.push_back(GetGray(0.25, 0.75, 0.75));
//	DisplaySpheres(spheres, Lines);
//	return 0;
//}
//int testCellList()
//{
//	Configuration a=ReadPos("3D/FCC");
//	for(;;)
//	{
//		double rc;
//		std::cin>>rc;
//		size_t count=0;
//		a.IterateThroughNeighbors(GeometryVector(a.GetDimension()), rc, [&count](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void {count++;});
//		std::cout<<count<<'\n';
//	}
//
//	return 0;
//}
/*
int DrawConfigs()
{
::AlsoWriteEPS=true;
Configuration c1=ReadPos("2D/Square");
c1.MoveParticle(c1.GetParticle(0), GeometryVector(0.5, 0.5));
::Plot("Square", MultiplicateStructure(c1, 3));
Configuration c2=ReadPos("2D/Paralelogram1.1");
c2.MoveParticle(c2.GetParticle(0), GeometryVector(0.5, 0.5));
::Plot("Paralelogram", MultiplicateStructure(c2, 3));
Configuration c3=ReadPos("2D/Triangle");
c3.MoveParticle(c3.GetParticle(0), GeometryVector(0.5, 0.5));
::Plot("Triangle", MultiplicateStructure(c3, 3));

return 0;
}
*/
//int TestMINOP()
//{
//	::Verbosity=3;
//	LoadStructure target("2D/Rectanglepi.pos");
//	GaussianCoreParameterization para(&target);
//	std::fstream parafile("tempPara.txt", std::fstream::in);
//	para.Read(parafile);
//	IsotropicPairPotential * pPot=para.GetPotential();
//	Configuration con=target.GetStructure();
//	FromStructure ts(con);
//	double Pressure=GetPressure(&ts, con.PeriodicVolume()/con.NumParticle(), *pPot);
//	RelaxStructure_GSL(con, *pPot, Pressure, 2, para.MinDistance());
//	std::cout<<"NLOPT~\n";
//	RelaxStructure_NLOPT(con, *pPot, Pressure, 2, para.MinDistance());
//	return 0;
//}

const double VirtualStructureSideLength=1.875;
const DimensionType VirtualStructureDimension=3;
class HCPSk : public LatticeSumSolver
{
public:
	HCPSk(void)
	{
		this->CurrentRc=0;
		this->OriginalDensity=1.0/std::pow(VirtualStructureSideLength, VirtualStructureDimension);
		this->NumParticle=1;
		this->Dimension=VirtualStructureDimension;
	}
	virtual const char * Tag(void)
	{
		return "HCP_Sk";
	}
	virtual ~HCPSk()
	{
	}
	virtual Configuration GetStructure(void)
	{
		Configuration result=GetUnitCubicBox(VirtualStructureDimension);
		result.Rescale(VirtualStructureSideLength);
		result.Insert("A", GeometryVector(VirtualStructureDimension));
		return result;
	}
	virtual void UpdateTerms(double NewRc)
	{
		this->Terms.clear();
		std::fstream ifile("HCPSkTerms.txt", std::fstream::in);
		for(;;)
		{
			LatticeTerm t;
			ifile>>t.distance;
			ifile>>t.number;
			if(ifile.eof()==true)
				break;
			if(t.distance>NewRc)
				break;
			this->Terms.push_back(t);
		}
		std::cout<<"Virtual lattice term updated, last term distance="<<this->Terms.back().distance<<'\n';
		this->CurrentRc=NewRc;
	}
};

#include <gsl/gsl_fit.h>
struct AuxType
{
	std::vector<Configuration> InitConfigs;
	ParameterSet * pParam;
	bool OutputStructureFactor;
};
double ObjectiveFunction(unsigned int n, const double *x, double *grad, void *aux)
{
	assert(grad==nullptr);
	AuxType * par =  reinterpret_cast<AuxType *> (aux);
	std::memcpy(par->pParam->Parameters, x, sizeof(double)*par->pParam->NumParameters());
	Potential * ppot = par->pParam->GetPotential();
	std::vector<Configuration> InherentConfigs;
	for(auto iter=par->InitConfigs.begin(); iter!=par->InitConfigs.end(); iter++)
	{
		Configuration c= *iter;
		RelaxStructure(c, *ppot, 1.0, par->pParam->MinDistance());
		InherentConfigs.push_back(c);
	}
	std::vector<GeometryVector> Sk;

	auto GetConfigsFunction = [&] (size_t i) ->Configuration
	{
		return InherentConfigs[i];
	};

	//supress outputing progress bar
	size_t TempVerbosity=1;
	std::swap(TempVerbosity, Verbosity);
	IsotropicStructureFactor(GetConfigsFunction, InherentConfigs.size(), 1.0, 1.0, Sk, 0.0);
	std::swap(TempVerbosity, Verbosity);
	if(par->OutputStructureFactor)
	{
		std::cout<<"Sk\n";
		for(auto iter=Sk.begin(); iter!=Sk.end(); iter++)
			std::cout<<iter->x[0]<<" \t"<<iter->x[1]<<'\n';
	}

	//std::vector<double> xx, y;
	//for(auto iter=Sk.begin(); iter!=Sk.end(); iter++)
	//{
	//	xx.push_back(iter->x[0]*iter->x[0]);
	//	y.push_back(iter->x[1]);
	//}
	//double c0, c1, cov00, cov01, cov11, sumsq;
	//gsl_fit_linear(&xx[0], 1, &y[0], 1, xx.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

	//delete ppot;
	//std::cout<<"p0="<<par->pParam->Parameters[0]<<", c0="<<c0<<", c1="<<c1<<", |residual|="<<std::sqrt(sumsq)<<'\n';
	//return c0;
	double sumError=1e-300;
	for(auto iter=Sk.begin(); iter!=Sk.end(); iter++)
	{
		double temp = iter->x[1] - 0.0025 * iter->x[0] * iter->x[0];

		if(iter==Sk.begin() || iter==Sk.end()-1 || iter==Sk.begin()+Sk.size()/2)
			std::cout<<temp<<" \t";

		sumError+=temp*temp;
	}
	double result=sumError/Sk.size();
	std::cout<<result<<'\n';
	return result;
}
#include "RandomSequentialAddition.h"
int InvStatMechDebug()
{	
	//DimensionType dim=2;
	//RandomGenerator gen(12345);
	//AuxType aux;
	//aux.OutputStructureFactor=false;
	//for(int i=0; i<5; i++)
	//{
	//	SpherePacking pk = GenerateRSAPacking(dim, 100, -3, 1000, 1000, 0.49, i, logfile);
	//	Configuration c(pk, "A");
	//	c.Resize(100.0);
	//	aux.InitConfigs.push_back(c);
	//}

	//struct OptimizationResult
	//{
	//	double MinF;
	//	ParameterSet * pParam;
	//};
	//std::vector<OptimizationResult> results;

	//for(int i=0; i<10; i++)
	//{
	//	std::cout<<"A new optimization:\n";
	//	aux.pParam = new GaussianCoreParameterization(dim, 0.3);
	//	aux.pParam->Evolve();

	//	ParameterSet & param=*aux.pParam;
	//	size_t dimension = param.NumParameters();
	//	double * variable = new double[dimension];
	//	//std::memcpy(variable, param.Parameters, param.NumParameters()*sizeof(double));
	//	double * LowerBounds = new double[dimension];
	//	double * UpperBounds = new double[dimension];
	//	double minf; // the minimum objective value, upon return
	//	nlopt_opt opt;
	//	opt = nlopt_create(NLOPT_LN_COBYLA, dimension); 
	//	nlopt_result status;
	//	status=nlopt_set_min_objective(opt, ObjectiveFunction, reinterpret_cast<void *>(&aux));
	//	param.GetBound(LowerBounds, UpperBounds);
	//	status=nlopt_set_lower_bounds(opt, LowerBounds);
	//	status=nlopt_set_upper_bounds(opt, UpperBounds);
	//	status=nlopt_set_ftol_abs(opt, 1e-15);
	//	for(int i=0; i<dimension; i++)
	//		variable[i]=LowerBounds[i]+gen.RandomDouble()*(UpperBounds[i]-LowerBounds[i]);
	//	status=nlopt_optimize(opt, variable, &minf);

	//	OptimizationResult temp;
	//	temp.MinF=minf;
	//	temp.pParam=aux.pParam;
	//	results.push_back(temp);


	//	delete [] variable, LowerBounds, UpperBounds;
	//}
	//auto CompareFunc = [] (const OptimizationResult & left, const OptimizationResult & right) ->bool
	//{
	//	return left.MinF<right.MinF;
	//};
	//std::sort(results.begin(), results.end(), CompareFunc);
	//aux.pParam=results[0].pParam;
	//std::cout<<"Lowest parameters:\n";
	//aux.pParam->Write(std::cout);
	//aux.OutputStructureFactor=true;
	//std::cout<<"testing objective fuction with minimum parameters\n";
	//ObjectiveFunction(aux.pParam->NumParameters(), aux.pParam->Parameters, nullptr, reinterpret_cast<void *>(&aux));
	//std::cout<<"testing objective fuction with more structures\n";
	//for(int i=0; i<50; i++)
	//{
	//	SpherePacking pk = GenerateRSAPacking(dim, 100, -3, 1000, 1000, 0.49, i, logfile);
	//	Configuration c(pk, "A");
	//	c.Resize(100.0);
	//	aux.InitConfigs.push_back(c);
	//}
	//ObjectiveFunction(aux.pParam->NumParameters(), aux.pParam->Parameters, nullptr, reinterpret_cast<void *>(&aux));



	return 0;
}
