#ifndef MONTECARLO_INCLUDED
#define MONTECARLO_INCLUDED

//General-Purpose Monte Carlo, not necessary for particles

#include "etc.h"
#include "RandomGenerator.h"
#include "Move.h"
#include <functional>


class CoolingSchedule
{
public:
	bool Continue;
	double Temperature;
	size_t NumTry;

	virtual void Report(double Energy) = 0;
	//when initialized, the variables are set to the 1st stage cooling.
	//after class SimulatedAnnealing performs these MC trials, it will call Report to report it's current energy
	//the class need to decide if another stage is needed, and update the variables
};

//General-Purpose Monte Carlo, not necessary for particles
template<typename System, typename AdditionalInfo> class MonteCarlo
{
protected:
	RandomGenerator gen;
public:

	double RecordEnergy;//this class keeps a log of the system energy, this is valid only if the potential is never changed
	size_t MoveCount;
	System * pSys;

	//Some trial moves, like particle displacement, has step size that needs to be adjusted.
	//The adjustment can be automatic. However, it must be turned off during sampling.
	//Set LockStepSize=true to turn off automatic adjustment of trial move parameters.
	bool LockStepSize;

	MonteCarlo(System & sys, int RandomSeed, double Energy) : gen(RandomSeed)
	{
		this->MoveCount=0;
		this->pSys= & sys;
		this->RecordEnergy=Energy;
		this->LockStepSize=false;
	}
	virtual ~MonteCarlo()
	{
	}

	virtual void Move(double Temperature, size_t Repeat, MCMove<System, AdditionalInfo> & move, AdditionalInfo add)
	{
		for(size_t i=0; i<Repeat; i++)
		{
			MoveCount++;
			double PreFactor=1.0;
			double dE=move.DeltaEnergy(*this->pSys, add, this->gen, PreFactor, LockStepSize, RecordEnergy);
			if(Temperature ==0)
			{
				if(dE>0)
					continue;
			}
			else if(gen.RandomDouble() > PreFactor*std::exp((-1)*dE/Temperature))
				continue;
				
			move.Accept(*this->pSys, add);
			this->RecordEnergy+=dE;
		}
	}
	virtual void Anneal(AdditionalInfo add, CoolingSchedule & cool, MCMove<System, AdditionalInfo> & move, std::function<void(double & Energy, double Temperature, System & sys)> CallBack)
	{
		while(cool.Continue)
		{
			this->Move(cool.Temperature, cool.NumTry, move, add);
			CallBack(this->RecordEnergy, cool.Temperature, *pSys);

			cool.Report(this->RecordEnergy);
		}
	}
};


template<typename System, typename AdditionalInfo, typename BinDecider> class WangLandauWorker;
//class for Wang-Landau Monte Carlo
//class BinDecider should have the following member functions:
//size_t NumBins()
//long GetBin(double E)
//double GetBinLowerBound(size_t NumBin)
//double GetBinUpperBound(size_t NumBin)
//The user is responsible for adjusting SIncrease
template<typename System, typename AdditionalInfo, typename BinDecider> class WangLandauMonteCarlo
{
	friend class WangLandauWorker <System, AdditionalInfo, BinDecider>;
private:
	RandomGenerator gen;
	BinDecider * pBinDecider;
	size_t NumBins;
public:
	std::vector<double> s;
	std::vector<size_t> Histogram;
	size_t MoveCount;

	bool LockStepSize;
	double RecordEnergy;
	System * pSys;
	double SIncrease;
	std::vector<System> vSys;
	WangLandauMonteCarlo(System & sys, int RandomSeed, double Energy, BinDecider & d, bool RecordConfigurationsPerBin = false) :
		gen(RandomSeed), s(d.NumBins(), 0.0), Histogram(d.NumBins(), 0), RecordEnergy(Energy)
	{
		this->MoveCount = 0;
		this->pSys = &sys;
		this->LockStepSize = false;
		this->pBinDecider = &d;
		this->NumBins = d.NumBins();
		if (d.GetBin(Energy) >= d.NumBins() || d.GetBin(Energy)<0)
		{
			std::cerr << "Error in WangLandauMonteCarlo : initial energy not within bounds! Bin=" << d.GetBin(Energy) << "\n";
			assert(false);
		}
		this->SIncrease=1.0;

		if (RecordConfigurationsPerBin)
			vSys.resize(d.NumBins());
	}
	void ClearHistogram(void)
	{
		for (size_t i = 0; i < NumBins; i++)
			Histogram[i] = 0;
	}
	std::vector<GeometryVector> GetHistogram(void)
	{
		std::vector<GeometryVector> result;
		for (int i = 0; i<NumBins; i++)
		{
			GeometryVector t(2);
			t.x[0] = (pBinDecider->GetBinLowerBound(i) + pBinDecider->GetBinUpperBound(i))*0.5;
			t.x[1] = this->Histogram[i];
			result.push_back(t);
		}
		return result;
	}
	std::vector<GeometryVector> GetEntropy(void)
	{
		std::vector<GeometryVector> result;
		for (int i = 0; i<NumBins; i++)
		{
			GeometryVector t(2);
			t.x[0] = (pBinDecider->GetBinLowerBound(i) + pBinDecider->GetBinUpperBound(i))*0.5;
			t.x[1] = this->s[i];
			result.push_back(t);
		}
		return result;
	}
	void Move(size_t Repeat, MCMove<System, AdditionalInfo> & move, AdditionalInfo add, std::function<double(const System & sys)> GetEnergyFunc)
	{
		for (size_t i = 0; i<Repeat; i++)
		{
			MoveCount++;
			double PreFactor = 1.0;
			double dE = move.DeltaEnergy(*this->pSys, add, this->gen, PreFactor, LockStepSize, RecordEnergy);
			long PrevBin = pBinDecider->GetBin(RecordEnergy);
			long AfterBin = pBinDecider->GetBin(RecordEnergy + dE);


			if ((PrevBin + 1 < 0) || (PrevBin > NumBins))
			{
				std::cerr << "Error in WangLandauMonteCarlo : PrevBin out of range!\n";
				std::cerr << "PrevBin=" << PrevBin << ", RecordEnergy=" << RecordEnergy << ", MoveCount=" << MoveCount << '\n';
				assert(false);
			}
			else if (PrevBin == -1)
				PrevBin = 0;
			else if (PrevBin == NumBins)
				PrevBin = NumBins - 1;
			if (AfterBin < 0 || AfterBin >= NumBins)
			{
				s[PrevBin] += SIncrease;
				Histogram[PrevBin]++;
			}
			else if (gen.RandomDouble() < PreFactor*std::exp(s[PrevBin] - s[AfterBin]))
			{
				move.Accept(*this->pSys, add);
				this->RecordEnergy+=dE;

				if (s[AfterBin] == 0 && vSys.size() != 0)
				{
					vSys[AfterBin] = (*this->pSys);
					//debug temp
					//ConfigurationPack pk("wlmc_sys");
					//pk.AddConfig(*this->pSys);
				}

				s[AfterBin]+=SIncrease;
				Histogram[AfterBin]++;
			}
			else
			{
				s[PrevBin] += SIncrease;
				Histogram[PrevBin]++;
			}
			if (MoveCount % 1000 == 0)
			{
				double newE = GetEnergyFunc(*pSys);
				if (newE<1e5 && std::abs((newE - RecordEnergy) / (RecordEnergy)) > 1e-8 && std::abs((newE - RecordEnergy)) > 1e-8)
					std::cerr << "Error in WangLandauMonteCarlo: Energy not consistent, newE=" << newE << ", record E=" << RecordEnergy << '\n';

				this->RecordEnergy = newE;
			}
		}
	}
	bool HistogramFlat(double tol)//return true if Histogram(E) for every possible E is not less than (1-tol)*(average)
	{
		double sum1 = 0.0, sumH = 0.0, minH = 1e300;
		for (size_t i = 0; i < NumBins; i++)
		{
			size_t temp = Histogram[i];
			if (temp > 0)
			{
				sum1 += 1.0;
				sumH += temp;
				minH = std::min(minH, (double)(temp));
			}
		}
		if (sum1 == 0.0)
			return false;
		else
			return minH >= (1.0 - tol)*sumH / sum1;
	}
};

//constructed from an WangLandauMonteCarlo object (object A), and uses A.Entropy and A.Histogram
//one can construct multiple WangLandauWorker from A, and then use them in parallel.
template<typename System, typename AdditionalInfo, typename BinDecider> class WangLandauWorker
{
private:
	double * pEntropy;
	size_t * pHistogram;
	size_t * pMoveCount;
	RandomGenerator mygen;
	BinDecider * pBinDecider;
	size_t NumBins;
public:
	bool LockStepSize;
	double RecordEnergy;
	System * pSys;
	double SIncrease;

	WangLandauWorker(System & sys, int RandomSeed, double Energy, WangLandauMonteCarlo<System, AdditionalInfo, BinDecider> & parent) :
		mygen(RandomSeed), pEntropy(parent.Entropy.data()), pHistogram(parent.Histogram.data()), RecordEnergy(Energy)
	{
		this->pMoveCount = &parent.MoveCount;
		this->pSys = &sys;
		this->LockStepSize = parent.LockStepSize;
		this->pBinDecider = parent.pBinDecider;
		this->NumBins = parent.NumBins;
		if (pBinDecider->GetBin(Energy) >= NumBins || pBinDecider->GetBin(Energy)<0)
		{
			std::cerr << "Error in WangLandauWorker : initial energy not within bounds! Bin=" << pBinDecider->GetBin(Energy) << "\n";
			assert(false);
		}
		this->SIncrease = parent.SIncrease;
	}

	void Move(size_t Repeat, MCMove<System, AdditionalInfo> & move, AdditionalInfo add, std::function<double(const System & sys)> GetEnergyFunc)
	{
		for (size_t i = 0; i<Repeat; i++)
		{
#pragma omp atomic
			(*pMoveCount)++;
			double PreFactor = 1.0;
			double dE = move.DeltaEnergy(*this->pSys, add, this->mygen, PreFactor, LockStepSize, RecordEnergy);
			long PrevBin = pBinDecider->GetBin(RecordEnergy);
			long AfterBin = pBinDecider->GetBin(RecordEnergy + dE);


			if ((PrevBin + 1 < 0) || (PrevBin > NumBins))
			{
				std::cerr << "Error in WangLandauWorker : PrevBin out of range!\n";
				std::cerr << "PrevBin=" << PrevBin << ", RecordEnergy=" << RecordEnergy << ", MoveCount=" << (*pMoveCount) << '\n';
				assert(false);
			}
			else if (PrevBin == -1)
				PrevBin = 0;
			else if (PrevBin == NumBins)
				PrevBin = NumBins - 1;
			if (AfterBin < 0 || AfterBin >= NumBins)
			{
#pragma omp atomic
				pEntropy[PrevBin] += SIncrease;
#pragma omp atomic
				pHistogram[PrevBin]++;
			}
			else if (mygen.RandomDouble() < PreFactor*std::exp(pEntropy[PrevBin] - pEntropy[AfterBin]))
			{
				move.Accept(*this->pSys, add);
				this->RecordEnergy += dE;
#pragma omp atomic
				pEntropy[AfterBin] += SIncrease;
#pragma omp atomic
				pHistogram[AfterBin]++;
			}
			else
			{
#pragma omp atomic
				pEntropy[PrevBin] += SIncrease;
#pragma omp atomic
				pHistogram[PrevBin]++;
			}
			if (*pMoveCount % 1000 == 0)
			{
				double newE = GetEnergyFunc(*pSys);
				if (newE<1e5 && std::abs((newE - RecordEnergy) / (RecordEnergy)) > 1e-8 && std::abs((newE - RecordEnergy)) > 1e-8)
					std::cerr << "Error in WangLandauMonteCarlo: Energy not consistent, newE=" << newE << ", record E=" << RecordEnergy << '\n';

				this->RecordEnergy = newE;
			}
		}
	}
};



const size_t IndividualMoveTime=25;
//const double SwapMovePropability=0.2;
template<typename System, typename AdditionalInfo> class ParallelTempering : public MonteCarlo<System, AdditionalInfo>
{
public:
	std::vector<double> vEnergyFactors;
	std::vector<MonteCarlo<System, AdditionalInfo> > vMonteCarlos;
	//one of the systems is given(others are allocated by this class), the pointer pSys points to it.
	size_t SwapSuccessCount, SwapFailCount;
	System * pOriginalSys;

	ParallelTempering<System, AdditionalInfo>(System & sys, int RandomSeed, size_t NumCopy, double Energy) : MonteCarlo<System, AdditionalInfo>(sys, RandomSeed, Energy)
	{
		if(NumCopy==0)
		{
			std::cerr<<"Error in Parallel Tempering : NumCopy is 0\n";
			assert(false);
		}
		{
			this->vEnergyFactors.clear();
			double nowRescale=1.0;
			for(size_t i=0; i<NumCopy; i++)
			{
				this->vEnergyFactors.push_back(nowRescale);
				nowRescale*=1.2;
			}
		}

		for(size_t i=0; i<this->vEnergyFactors.size(); i++)
		{
			if(i==0)
				this->vMonteCarlos.push_back(MonteCarlo<System, AdditionalInfo>(sys, (RandomSeed+i), Energy ));
			else
				this->vMonteCarlos.push_back(MonteCarlo<System, AdditionalInfo>(*(new System(sys)), (RandomSeed+i), Energy ));
			this->vMonteCarlos[i].RecordEnergy=Energy;
		}
		this->pSys = & sys;
		this->pOriginalSys= & sys;
		this->SwapFailCount=0;
		this->SwapSuccessCount=0;
	}
	~ParallelTempering()
	{
		//debug temp
		//std::cerr<<"Parallel Tempering destructed, Fail times:"<<this->SwapFailCount<<", Success Count:"<<this->SwapSuccessCount<<'\n';

		for(size_t i=0; i<this->vMonteCarlos.size(); i++)
			if(this->vMonteCarlos[i].pSys!= this->pOriginalSys)
				delete this->vMonteCarlos[i].pSys;
	}

	void IndividualMove(double Temperature, size_t Repeat, std::vector<MCMove<System, AdditionalInfo> * > & moves, AdditionalInfo add)
	{
		size_t stop = this->vEnergyFactors.size();
		assert(this->vMonteCarlos.size() == stop);
		assert(moves.size() == stop);

		size_t NumThreads=PTParallelConfigurationsPerSystem;
		if(stop < NumThreads)
			NumThreads=stop;

		#pragma omp parallel for num_threads(NumThreads)
		for(signed long i=0; i<stop; i++)
			this->vMonteCarlos[i].Move(Temperature/this->vEnergyFactors[i], Repeat, *moves[i], add);
	}
	virtual void Move(double Temperature, size_t Repeat, std::vector<MCMove<System, AdditionalInfo> * > & moves, AdditionalInfo add)
	{
		//sync this->RecordEnergy
		for(auto iter = this->vMonteCarlos.begin(); iter!=this->vMonteCarlos.end(); iter++)
		{
			if(iter->pSys == this->pSys)
			{
				iter->RecordEnergy=this->RecordEnergy;
				break;
			}
		}

		for(size_t i=0; i<Repeat; i++)
		{
			if(this->gen.RandomDouble() < 0.9)
			{
				//swap moves
				for(size_t i=0; i<this->vMonteCarlos.size(); i++)
				{
					size_t m=static_cast<size_t>(this->gen.RandomDouble()*100000000)%(this->vMonteCarlos.size()-1);
					size_t n=m+1;
					if(n>=this->vMonteCarlos.size())
					{
						std::cerr<<"error in ParallelTempering : index out of range\n";
						continue;
					}

					double DeltaE=this->vMonteCarlos[n].RecordEnergy*this->vEnergyFactors[m]+this->vMonteCarlos[m].RecordEnergy*this->vEnergyFactors[n]-this->vMonteCarlos[n].RecordEnergy*this->vEnergyFactors[n]-this->vMonteCarlos[m].RecordEnergy*this->vEnergyFactors[m];

					if(Temperature ==0)
					{
						if(DeltaE>0)
						{
							this->SwapFailCount++;
							continue;
						}
					}
					else if(this->gen.RandomDouble() > std::exp((-1)*DeltaE/Temperature))
					{
						this->SwapFailCount++;
						continue;
					}

					std::swap(this->vMonteCarlos[m].pSys, this->vMonteCarlos[n].pSys);
					std::swap(this->vMonteCarlos[m].RecordEnergy, this->vMonteCarlos[n].RecordEnergy);
					this->SwapSuccessCount++;

					//debug temp
					if(this->pSys!=this->vMonteCarlos.back().pSys)
					{
						this->pSys=this->vMonteCarlos.back().pSys;
						//std::cerr<<"pSys changed!";
					}

				}

				//modify energy factors if necessary
				size_t count=this->SwapFailCount+this->SwapSuccessCount;
				if(count>1000 && this->vEnergyFactors.size()>1)
				{
					//debug temp
					//std::cout<<"Swap success:"<<this->SwapSuccessCount<<", fail:"<<this->SwapFailCount<<'\n';
					if(static_cast<double>(this->SwapSuccessCount)/count<0.1)
					{
						double factor=std::pow(this->vEnergyFactors[1], 0.7);
						for(size_t i=0; i<this->vEnergyFactors.size(); i++)
							this->vEnergyFactors[i]=std::pow(factor, static_cast<double>(i));

						//debug temp
						std::cout<<"Energy factor has been changed to:"<<factor<<'\n';
					}
					this->SwapFailCount=0;
					this->SwapSuccessCount=0;
				}
			}
			else
			{
				this->IndividualMove(Temperature, IndividualMoveTime, moves, add);
				i+=IndividualMoveTime;
			}
		}
		this->MoveCount += Repeat;

		//sync this->RecordEnergy
		for(auto iter = this->vMonteCarlos.begin(); iter!=this->vMonteCarlos.end(); iter++)
		{
			if(iter->pSys == this->pSys)
			{
				this->RecordEnergy=iter->RecordEnergy;
				break;
			}
		}
	}
	virtual void Move(double Temperature, size_t Repeat, MCMove<System, AdditionalInfo>  & move, AdditionalInfo add)
	{
		std::vector<MCMove<System, AdditionalInfo> * > moves;
		for(size_t i=0; i<this->vMonteCarlos.size(); i++)
			moves.push_back(move.clone());

		this->Move(Temperature, Repeat, moves, add);
		//delete stuff
		for(auto iter=moves.begin(); iter!=moves.end(); iter++)
			delete *iter;
	}
	virtual void Anneal(AdditionalInfo add, CoolingSchedule & cool, MCMove<System, AdditionalInfo> & move, std::function<void(double & Energy, double Temperature, System & sys)> CallBack)
	{
		std::vector<MCMove<System, AdditionalInfo> * > moves;
		for(size_t i=0; i<this->vMonteCarlos.size(); i++)
			moves.push_back(move.clone());

		while(cool.Continue)
		{
			this->Move(cool.Temperature, cool.NumTry, moves, add);
			CallBack(this->RecordEnergy, cool.Temperature, *this->pSys);

			//report the average energy of all configurations to the cooling schedule
			double AverageEnergy=0;
			for(size_t i=0; i<this->vMonteCarlos.size(); i++)
				AverageEnergy+=this->vEnergyFactors[i]*this->vMonteCarlos[i].RecordEnergy;
			AverageEnergy/=this->vMonteCarlos.size();
			cool.Report(AverageEnergy);

			//debug temp
			//for(size_t i=0; i<cool.NumTry; i++)
			//{
			//	this->Move(cool.Temperature, 1, move, add);
			//	CallBack(this->RecordEnergy, cool.Temperature, *pSys);
			//}
			//cool.Report(this->RecordEnergy);
		}
		//delete stuff
		for(auto iter=moves.begin(); iter!=moves.end(); iter++)
			delete *iter;

		//the one with lowest temperature is returned
		if(this->pOriginalSys!=this->pSys)
			(*this->pOriginalSys)=(*this->vMonteCarlos.back().pSys);
	}
};


const size_t StageSize=1000;
class ExpCool : public CoolingSchedule
{
	double StageAlpha;
	double Tcut;
	signed long ZeroMove;
public:
	ExpCool(double Alpha, double Tinit, double Tcut, long ZeroMove) : Tcut(Tcut), ZeroMove(ZeroMove)
	{
		this->StageAlpha=std::pow(Alpha, static_cast<double>(StageSize));
		this->Continue=true;
		this->Temperature=Tinit;
		this->NumTry=StageSize;
	}
	virtual void Report(double Energy)
	{
		//////////////////////////
		if(this->Temperature==0)
		{
			this->ZeroMove-=this->NumTry;
			if(this->ZeroMove<=0)
				this->Continue=false;
		}
		else
		{
			this->Temperature*=this->StageAlpha;
			if(this->Temperature<this->Tcut)
			{
				this->Temperature=0;
				if(this->NumTry>this->ZeroMove)
					this->NumTry=this->ZeroMove;
			}
		}
	}
};



const size_t ReportTime=100;
class ZeroCool : public CoolingSchedule
{
	size_t NumReport;

public:
	ZeroCool(long ZeroMoveTime)
	{
		this->Temperature=0.0;
		this->NumTry=ZeroMoveTime/ReportTime;
		this->Continue=true;
		this->NumReport=0;
	}
	virtual void Report(double Energy)
	{
		this->NumReport++;
		if(this->NumReport>ReportTime)
			this->Continue=false;
	}
};


#include "AnalyzeData.h"
const size_t DefaultSampleInterval = 100;
const double DefaultVsCoeff=500;
const size_t DefaultInitAnalyzeSize=200;
class ThermoCool : public CoolingSchedule
{
	size_t SampleInterval;
	double VsCoeff;
	size_t InitAnalyzeSize;

	double Vs;
	double Tcut;
	signed long ZeroMove;//deprecated parameter

	size_t AnalyzeSize;

	std::vector<double> data;
public:
	ThermoCool(double Alpha, double Tinit, double Tcut, long ZeroMove, 
		size_t SampleInterval=DefaultSampleInterval, double VsCoeff=DefaultVsCoeff, 
		size_t InitAnalyzeSize=DefaultInitAnalyzeSize) : Tcut(Tcut), ZeroMove(ZeroMove), 
		SampleInterval(SampleInterval), VsCoeff(VsCoeff), InitAnalyzeSize(InitAnalyzeSize)
	{
		this->Vs=(1.0-Alpha)*SampleInterval*VsCoeff;
		this->Continue=true;
		this->Temperature=Tinit;
		this->NumTry=SampleInterval;

		this->AnalyzeSize=InitAnalyzeSize;
	}
	virtual void Report(double Energy)
	{
		//////////////////////////
		if(this->Temperature==0)
		{
			this->Continue=false;
		}
		else
		{
			data.push_back(Energy);

			if(data.size()>=AnalyzeSize)
			{
				double Average, error, Corrtime;
				AnalyzeData(data, Average, error, Corrtime);
				double H2c=error*error*data.size()/(1+2*Corrtime);//this is T^2*C_v
				double dlnT=this->Vs*this->Temperature/Corrtime/std::sqrt(H2c);

					//debug temp
					//std::fstream file("ErrorLog.txt", std::fstream::out);
					//for(auto iter=data.begin(); iter!=data.end(); iter++)
					//	file<<(*iter)<<'\n';
					//file.close();
					//exit(0);
				if( (!(dlnT>0)) )
				{
					if(this->AnalyzeSize < 5*InitAnalyzeSize)
					{
						this->Vs*=2;
						this->AnalyzeSize*=2;
					}
					this->Temperature/=std::exp(0.000001);//cool very cautiously when can't analyze data
				}
				else if(dlnT>0.1)
				{
					//std::cout<<"Warning in ThermoCool : Temperature dropping too quickly, decrease sample size.\n";
					if(this->AnalyzeSize>InitAnalyzeSize/2 && this->AnalyzeSize>1)
					{
						this->Vs/=2;
						this->AnalyzeSize/=2;
					}
					this->Temperature/=std::exp(0.1);
				}
				else
				{
					this->Temperature/=std::exp(dlnT);
				}

				data.clear();
			}

			if(this->Temperature<this->Tcut)
			{
				this->Temperature=0;
				if(this->ZeroMove==0)
					this->Continue=false;
				else
					this->NumTry=this->ZeroMove;
			}
		}
	}
};

#endif
