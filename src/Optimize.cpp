#include "Optimize.h"
#include "GeometryVector.h"
#include "Potential.h"
#include "Plots.h"
#include "Solvers.h"
#include "MC_System.h"
#include "StructureOptimization.h"

#include <omp.h>

#include <nlopt.hpp>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cmath>

//#define logfile (std::cerr)
//#define RecordCompetitors

struct AuxType//struct used for optimization
{
	const std::vector<LatticeSumSolver *> * pvpSolvers;
	ParameterSet * pParam;
	double TargetStructureVolume;
	size_t NumEvaluation;
	bool AlsoConsiderElastic;
};
void AllDerivatives(double *result, Potential * pPot, Configuration & tar, double Pressure)
{
	DimensionType dim=tar.GetDimension();
	size_t NumP=tar.NumParticle();

	pPot->SetConfiguration(tar);
	for(size_t i=0; i<NumP; i++)
	{
		//Configuration::particle * pA=tar.GetParticle(i);
		GeometryVector Force;
		pPot->Force(Force, i);

		//derivatives of atom coordinates
		for(DimensionType j=0; j<dim; j++)
		{
			result[dim*dim+i*dim+j] = (-1.0)*Force.Dot(tar.GetBasisVector(j));
		}
	}
	pPot->EnergyDerivativeToBasisVectors(result, Pressure);
	return;
}


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
double get_det(gsl_matrix * mat_ptr); //defined in GeometryVector.cpp
double LowestEigenValueOfStiffness(Configuration list, Potential & pot, double Pressure)
{
	//RelaxStructure(list, pot, Pressure, 0.0);
	DimensionType dim=list.GetDimension();
	size_t dimTensor=dim*(dim+1)/2;//dimension of stiffness tensor

	gsl_matrix * Stiffness = gsl_matrix_calloc(dimTensor, dimTensor);
	size_t Count1=0, Count2=0;//indexes of stiffness tensor
	for(DimensionType i=0; i<dim; i++)
	{
		for(DimensionType j=i; j<dim; j++)
		{
			for(DimensionType m=0; m<dim; m++)
			{
				for(DimensionType n=m; n<dim; n++)
				{
					if(Count1>=Count2)//this element not calculated yet
					{
						double temp=ElasticConstant(pot, list, i, j, m, n, Pressure);
						gsl_matrix_set(Stiffness, Count1, Count2, temp);
						gsl_matrix_set(Stiffness, Count2, Count1, temp);
					}
					Count2++;
				}
			}
			assert(Count2==dimTensor);
			Count2=0;
			Count1++;
		}
	}
	assert(Count1==dimTensor);
	gsl_eigen_symm_workspace * wor = gsl_eigen_symm_alloc(dimTensor);
	gsl_vector * eig = gsl_vector_calloc(dimTensor);

	gsl_eigen_symm(Stiffness, eig, wor);
	double LowestEigen = ::MaxEnergy;
	for(size_t i=0; i<dimTensor; i++)
	{
		double temp=gsl_vector_get(eig, i);
		if(temp<LowestEigen)
			LowestEigen=temp;
	}

	gsl_vector_free(eig);
	gsl_eigen_symm_free(wor);
	gsl_matrix_free(Stiffness);

	return LowestEigen;
}
double LowestEigenValueOfDeformations(Configuration list, Potential & pot, double Pressure)
{
	const double TrialSize=1e-4;
	//RelaxStructure(list, pot, Pressure, 0.0);
	DimensionType dim=list.GetDimension();
	size_t nbr=list.NumParticle();
	size_t dimVector=dim*(dim+nbr);
	double * temp1= new double [dimVector];
	double * temp2= new double [dimVector];
	size_t dimTensor=(dim*(dim+1))/2+(nbr-1)*dim;
	std::vector<DimensionType> ms, ns;
	for(DimensionType m=0; m<dim; m++)
	{
		for(DimensionType n=m; n<dim; n++)
		{
			ms.push_back(m);
			ns.push_back(n);
		}
	}
	assert(ms.size()==(dim*(dim+1))/2);

	gsl_matrix * Stiffness = gsl_matrix_calloc(dimTensor, dimTensor);
	AllDerivatives(temp1, &pot, list, Pressure);
	GeometryVector bas[ ::MaxDimension];
	for(DimensionType i=0; i<dim; i++)
		bas[i]=list.GetBasisVector(i);

	for(size_t i=0; i<ms.size(); i++)
	{
		DimensionType m=ms[i], n=ns[i];
		bas[n].x[m]+=TrialSize;
		list.ChangeBasisVector(bas);
		AllDerivatives(temp2, &pot, list, Pressure);
		//restore
		bas[n].x[m]-=TrialSize;
		list.ChangeBasisVector(bas);

		for(size_t j=0; j<ms.size(); j++)
		{
			size_t index=ns[j]*dim+ms[j];
			gsl_matrix_set(Stiffness, i, j, (temp2[index]-temp1[index])/TrialSize );
		}
		for(size_t j=0; j<dim*(nbr-1); j++)
		{
			gsl_matrix_set(Stiffness, i, ms.size()+j, (temp2[dim*dim+j]-temp1[dim*dim+j])/TrialSize );
			//the last particle has a trivial translation degree of freedom
		}
	}
	for(size_t n=0; n<nbr-1; n++)
	{
		//Configuration::particle * pA=list.GetParticle(n);
		GeometryVector aRel=list.GetRelativeCoordinates(n);

		for(DimensionType m=0; m<dim; m++)
		{
			aRel.x[m]+=TrialSize;
			list.MoveParticle(n, aRel);
			AllDerivatives(temp2, &pot, list, Pressure);
			//restore
			aRel.x[m]-=TrialSize;
			list.MoveParticle(n, aRel);

			for(size_t j=0; j<ms.size(); j++)
			{
				size_t index=ns[j]*dim+ms[j];
				gsl_matrix_set(Stiffness, ms.size()+n*dim+m, j, (temp2[index]-temp1[index])/TrialSize );
			}
			for(size_t j=0; j<dim*(nbr-1); j++)
			{
				gsl_matrix_set(Stiffness, ms.size()+n*dim+m, ms.size()+j, (temp2[dim*dim+j]-temp1[dim*dim+j])/TrialSize );
				//the last particle has a trivial translation degree of freedom
			}
		}
	}


	gsl_eigen_symm_workspace * wor = gsl_eigen_symm_alloc(dimTensor);
	gsl_vector * eig = gsl_vector_calloc(dimTensor);

	gsl_eigen_symm(Stiffness, eig, wor);

	std::vector<double> eigens;
	for(size_t i=0; i<dimTensor; i++)
	{
		double temp=gsl_vector_get(eig, i);
		eigens.push_back(temp);
		//debug temp
		//std::cout<<"Eig"<<temp<<'\n';
	}
	std::nth_element(eigens.begin(), eigens.begin(), eigens.end());
	double LowestEigen=eigens[0];

	gsl_vector_free(eig);
	gsl_eigen_symm_free(wor);
	gsl_matrix_free(Stiffness);
	delete [] temp1;
	delete [] temp2;

	return LowestEigen;
}
void DeformationsMatrixPrincipleMinorDeterminants(Configuration list, Potential & pot, double Pressure, double * results)
{
	const double TrialSize=1e-4;
	//RelaxStructure(list, pot, Pressure, 0.0);
	DimensionType dim=list.GetDimension();
	size_t nbr=list.NumParticle();
	size_t dimVector=dim*(dim+nbr);
	double * temp1= new double [dimVector];
	double * temp2= new double [dimVector];
	size_t dimTensor=(dim*(dim+1))/2+(nbr-1)*dim;
	std::vector<DimensionType> ms, ns;
	for(DimensionType m=0; m<dim; m++)
	{
		for(DimensionType n=m; n<dim; n++)
		{
			ms.push_back(m);
			ns.push_back(n);
		}
	}
	assert(ms.size()==(dim*(dim+1))/2);

	gsl_matrix * Stiffness = gsl_matrix_calloc(dimTensor, dimTensor);
	AllDerivatives(temp1, &pot, list, Pressure);
	GeometryVector bas[ ::MaxDimension];
	for(DimensionType i=0; i<dim; i++)
		bas[i]=list.GetBasisVector(i);

	for(size_t i=0; i<ms.size(); i++)
	{
		DimensionType m=ms[i], n=ns[i];
		bas[n].x[m]+=TrialSize;
		list.ChangeBasisVector(bas);
		AllDerivatives(temp2, &pot, list, Pressure);
		//restore
		bas[n].x[m]-=TrialSize;
		list.ChangeBasisVector(bas);

		for(size_t j=0; j<ms.size(); j++)
		{
			size_t index=ns[j]*dim+ms[j];
			gsl_matrix_set(Stiffness, i, j, (temp2[index]-temp1[index])/TrialSize );
		}
		for(size_t j=0; j<dim*(nbr-1); j++)
		{
			gsl_matrix_set(Stiffness, i, ms.size()+j, (temp2[dim*dim+j]-temp1[dim*dim+j])/TrialSize );
			//the last particle has a trivial translation degree of freedom
		}
	}
	for(size_t n=0; n<nbr-1; n++)
	{
		//Configuration::particle * pA=list.GetParticle(n);
		GeometryVector aRel=list.GetRelativeCoordinates(n);

		for(DimensionType m=0; m<dim; m++)
		{
			aRel.x[m]+=TrialSize;
			list.MoveParticle(n, aRel);
			AllDerivatives(temp2, &pot, list, Pressure);
			//restore
			aRel.x[m]-=TrialSize;
			list.MoveParticle(n, aRel);

			for(size_t j=0; j<ms.size(); j++)
			{
				size_t index=ns[j]*dim+ms[j];
				gsl_matrix_set(Stiffness, ms.size()+n*dim+m, j, (temp2[index]-temp1[index])/TrialSize );
			}
			for(size_t j=0; j<dim*(nbr-1); j++)
			{
				gsl_matrix_set(Stiffness, ms.size()+n*dim+m, ms.size()+j, (temp2[dim*dim+j]-temp1[dim*dim+j])/TrialSize );
				//the last particle has a trivial translation degree of freedom
			}
		}
	}

	results[0]=gsl_matrix_get(Stiffness, 0, 0);
	for(size_t i=1; i<dimTensor; i++)
	{
		gsl_matrix * minor = gsl_matrix_calloc(i+1, i+1);
		for(size_t j=0; j<i+1; j++)
			for(size_t k=0; k<i+1; k++)
				gsl_matrix_set(minor, j, k, gsl_matrix_get(Stiffness, j, k));
		results[i] = get_det(minor);
		gsl_matrix_free(minor);
	}

	gsl_matrix_free(Stiffness);
	delete [] temp1;
	delete [] temp2;

	logfile<<"Deformation tensor minor determinants are:";
	for(size_t i=0; i<dimTensor; i++)
	{
		logfile<<results[i]<<" \t";
	}
	logfile<<'\n';

	return;
}

double NearestNeighborDistance2Volume(LatticeSumSolver * Solver, double NearestNeighborDistance)
{
	double R=1.0;
	while(Solver->Terms.size()<2)
	{
		Solver->UpdateTerms(R);
		R*=2;
	}
	double ReScale = std::pow(NearestNeighborDistance/Solver->Terms[1].distance, static_cast<double>(Solver->Dimension));
	double result = (1.0/Solver->OriginalDensity)*ReScale;
	return result;
}

const double DeltaVolumeAddOne = 1.0 + 1e-9;
double GetPressure(LatticeSumSolver * l, double volume, IsotropicPairPotential & pot)//Get the equilibrium pressure of this structure at this specific volume
{
	double v1=volume*DeltaVolumeAddOne;
	double v2=volume/DeltaVolumeAddOne;
	double e1=l->LatticeSum(v1, pot);
	double e2=l->LatticeSum(v2, pot);
	double result=(e1-e2)/(v2-v1);
	if(result<0)
		result=0;
	return result;
}

//make sure that every LatticeSumSoler in that vector doesn't need to go to MinDistance, decrease the pressure if necessary
void CheckPressure(double & Pressure, IsotropicPairPotential & pot, const std::vector<LatticeSumSolver *> & vpSolvers, double MinDistance)
{
	for(auto iter = vpSolvers.begin(); iter!=vpSolvers.end(); iter++)
	{
		double temp = GetPressure(*iter, NearestNeighborDistance2Volume(*iter, MinDistance), pot);
		if(temp<Pressure)
		{
			//std::cerr<<"In GetPressure, "<<(*iter)->Tag()<<" would reach MinDistance if Pressure="<<temp<<" Pressure is:"<<Pressure<<'\n';
			Pressure=temp;
		}
	}
}
double GetPressure(LatticeSumSolver * l, double volume, IsotropicPairPotential & pot, const std::vector<LatticeSumSolver *> & vpSolvers, double MinDistance)//with additional constraint that no nearest neighbor of any structure is closer than MinDistance
{
	double Pressure = GetPressure(l, volume, pot);
	CheckPressure(Pressure, pot, vpSolvers, MinDistance);
	return Pressure;
}



namespace EnthalpyDifference
{
	double EnthalpyDifference(ParameterSet & param, const std::vector<LatticeSumSolver *> & Solver, double TargetStructureVolume, double Pressure, Configuration * pLowestCompetitor)
	{
		size_t Threads = ::OptimizationNumThreadsPerTrial;
		if(Threads == 0)
			Threads = omp_get_num_procs();
		if(Threads > (Solver.size()-1) )
			Threads = Solver.size()-1;

		std::vector<double> dEnthalpies;
		dEnthalpies.resize(Solver.size());

		IsotropicPairPotential * pPot = param.GetPotential();
		double RawSearchStepSize = param.SharpestRelativeSize()/pPot->Dimension+1;

		double MinDistance = param.MinDistance();

		//calculation
		double globalLowest = ::MaxEnergy/2;//a big number
		double globalLowestVolume;
		size_t globalnLowestSolver;

		/*
		Configuration Target = Solver[0]->GetStructure();
		double lTarget=(pPrecise.Energy(Target) + Pressure*Target.PeriodicVolume())/Target.NumParticle();
		*/

		Configuration RelaxTarget = Solver[0]->GetStructure();
		RelaxStructure(RelaxTarget, *pPot, Pressure, MinDistance);
		double lTarget=(pPot->Energy(RelaxTarget) + Pressure*RelaxTarget.PeriodicVolume())/RelaxTarget.NumParticle();
		FromStructure RelaxTargetSolver(RelaxTarget);
		RelaxTargetSolver.UpdateTerms(5.0*std::pow(RelaxTarget.PeriodicVolume()/RelaxTarget.NumParticle(), 1.0/RelaxTarget.GetDimension()));

		omp_lock_t lock;
		omp_init_lock(&lock);
		omp_lock_t lock2;
		omp_init_lock(&lock2);
		size_t stop=Solver.size();
//#pragma	omp parallel default(shared) num_threads(Threads)
		{
			double threadLowest = ::MaxEnergy/2;//a big number
			double threadLowestVolume=100.0;
			size_t threadnLowestSolver=0;
		IsotropicPairPotential * ThreadpPot = param.GetPotential();
//#pragma omp for nowait schedule(dynamic)
			for(signed long i=0; i<stop; i++)
			{
				Configuration comp=Solver[i]->GetStructure();
				try
				{
					RelaxStructure(comp, *ThreadpPot, Pressure, MinDistance);
				}
				catch(NotANumberFound & a)
				{
					std::cerr<<"output parameterization:\n";
					param.Write(std::cerr);
					std::cerr<<"Pressure="<<Pressure<<'\n';
					throw;
				}

				FromStructure temp(comp);
				////////////////////////////
				if(SameSolver(&RelaxTargetSolver, &temp)==false)
				{
					double dEnthalpy=(ThreadpPot->Energy(comp) + Pressure*comp.PeriodicVolume())/comp.NumParticle() - lTarget;
					dEnthalpies[i]=dEnthalpy;
					if(dEnthalpy < threadLowest)
					{
						threadLowest=dEnthalpy;
						threadLowestVolume=comp.PeriodicVolume()/comp.NumParticle();
						threadnLowestSolver=i;
					}
					//debug temp
					/*
					omp_set_lock(&lock2);
					{
					logfile<<"Competitor"<<i<<"Give delta Energy of"<<dEnthalpy<<'\n';
					}
					omp_unset_lock(&lock2);
					*/
				}
			}

			omp_set_lock(&lock);
			{
				if(threadLowest<globalLowest)
				{
					globalLowest=threadLowest;
					globalLowestVolume=threadLowestVolume;
					globalnLowestSolver=threadnLowestSolver;
				}
			}
			omp_unset_lock(&lock);
			delete ThreadpPot;
		}

		if(pLowestCompetitor != nullptr)
		{
			(*pLowestCompetitor)=Solver[globalnLowestSolver]->GetStructure();
			try
			{
				RelaxStructure((*pLowestCompetitor), *pPot, Pressure, MinDistance);
			}
			catch(NotANumberFound & a)
			{
				std::cerr<<"output parameterization:\n";
				param.Write(std::cerr);
				throw;
			}
		}


		omp_destroy_lock(&lock);
		omp_destroy_lock(&lock2);
		delete pPot;
		logfile<<"dEnthalpies are:";
		for(auto iter=dEnthalpies.begin(); iter!=dEnthalpies.end(); iter++)
			logfile<<*iter<<" \t";
		logfile<<'\n';

		return (-1)*globalLowest;
	}

	//no force on target
	void Equality_MultiConstraint(unsigned m, double *result, unsigned n, const double* x, double* grad, void* aux)
	{
		assert(grad==nullptr);
		AuxType * par =  reinterpret_cast<AuxType *> (aux);
		std::memcpy(par->pParam->Parameters, x, sizeof(double)*par->pParam->NumParameters());
		IsotropicPairPotential * pPot = par->pParam->GetPotential();
		Configuration tar = par->pvpSolvers->at(0)->GetStructure();
		double Pressure = GetPressure(par->pvpSolvers->at(0), par->TargetStructureVolume, *pPot, *(par->pvpSolvers), par->pParam->MinDistance() );

		AllDerivatives(result, pPot, tar, Pressure);
		delete pPot;

		logfile<<"Equality Constraints are:";
		for(size_t i=0; i<m; i++)
			logfile<<result[i]<<" \t";
		logfile<<'\n';
	}

	bool LocallyStable(Configuration tar, Potential & Pot, double Pressure, double MinDistance)
	{
		FromStructure t1(tar);
		RelaxStructure(tar, Pot, Pressure, MinDistance);
		FromStructure t2(tar);

		//debug temp
		//std::cout<<"in LOcallyStable, d="<<SolverDistance(&t1, &t2)<<'\n';

		return SameSolver(&t1, &t2);
	}

	//distance between target and local minimum is 0
	double EqualityConstraint(unsigned int n, const double *x, double *grad, void *aux)
	{
		for (size_t i = 0; i<n; i++)
		{
			if (x[i] != x[i])
			{
				logfile << "Error in EnthalpyDifference::EqualityConstraint : NaN detected\n";
				return 0;
			}
		}
		if (x[0] < 0)
		{
			std::cerr << "Error in EnthalpyDifference::EqualityConstraint : x[0]=" << x[0] << '\n';
			return 0;
		}
		assert(grad == nullptr);
		AuxType * par =  reinterpret_cast<AuxType *> (aux);
		std::memcpy(par->pParam->Parameters, x, sizeof(double)*par->pParam->NumParameters());
		IsotropicPairPotential * pPot = par->pParam->GetPotential();
		Configuration tar = par->pvpSolvers->at(0)->GetStructure();
		double Pressure = GetPressure(par->pvpSolvers->at(0), par->TargetStructureVolume, *pPot, *(par->pvpSolvers), par->pParam->MinDistance() );
		double tarPressure = GetPressure(par->pvpSolvers->at(0), par->TargetStructureVolume, *pPot);
		if(tarPressure!=Pressure)
		{
			delete pPot;
			return ::MaxEnergy;
		}
		double PrevSpecificVolume = tar.PeriodicVolume() / tar.NumParticle();
		RelaxStructure(tar, *pPot, Pressure, par->pParam->MinDistance());
		double AfterSpecificVolume = tar.PeriodicVolume() / tar.NumParticle();
		double dtemp = (PrevSpecificVolume - AfterSpecificVolume) / (PrevSpecificVolume + AfterSpecificVolume);
		FromStructure temp(tar);
		double result=SolverDistance(&temp, par->pvpSolvers->at(0))+dtemp*dtemp;
		par->NumEvaluation++;
		logfile<<"at evaluation"<<par->NumEvaluation<<", Constraint returns:"<<result<<'\n';
		delete pPot;
		return result;
	}
	// positive Equality_MultiConstraint, plus
	// negative Equality_MultiConstraint, plus
	// negative Hessian matrix principle minor determinants
	void InEquality_MultiConstraint(unsigned m, double *result, unsigned n, const double* x, double* grad, void* aux)
	{
		for(size_t i=0; i<n; i++)
		{
			if(x[i]!=x[i])
			{
				logfile<<"Error in EnthalpyDifference::ObjectiveFunction : NaN detected\n";
				for(size_t i=0; i<m; i++)
					result[i]=0.0;
				return;
			}
		}
		if (x[0] < 0)
		{
			std::cerr << "Error in EnthalpyDifference::InEquality_MultiConstraint : x[0]=" << x[0] << '\n';
			for(size_t i=0; i<m; i++)
				result[i]=0.0;
			return;
		}
		AuxType * par = reinterpret_cast<AuxType *> (aux);
		std::memcpy(par->pParam->Parameters, x, sizeof(double)*par->pParam->NumParameters());
		IsotropicPairPotential * pPot = par->pParam->GetPotential();
		Configuration tar = par->pvpSolvers->at(0)->GetStructure();
		double Pressure = GetPressure(par->pvpSolvers->at(0), par->TargetStructureVolume, *pPot, *(par->pvpSolvers), par->pParam->MinDistance() );
		size_t NumEqualityConstraint=tar.GetDimension()*(tar.NumParticle()+tar.GetDimension());
		Equality_MultiConstraint(NumEqualityConstraint, result, n, x, grad, aux);
		for(size_t i=0; i<NumEqualityConstraint; i++)
			result[NumEqualityConstraint+i]=(-1.0)*result[i];

		DeformationsMatrixPrincipleMinorDeterminants(tar, *pPot, Pressure, result+2*NumEqualityConstraint);
		delete pPot;
		//allow some error in forces
		for(size_t i=0; i<m; i++)
			result[i]-=2e-3;
		//require positive minor determinants
		for(double * p=result+2*NumEqualityConstraint; p<result+m; p++)
			(*p)=(*p)*(-1.0)+0.01;

		logfile<<"InEquality Constraints are:";
		for(size_t i=0; i<m; i++)
			logfile<<result[i]<<" ";
		logfile<<'\n';
	}
	double ObjectiveFunction(unsigned int n, const double *x, double *grad, void *aux)//Enthalpy Difference at certail Pressure
	{
		for(size_t i=0; i<n; i++)
		{
			if(x[i]!=x[i])
			{
				logfile<<"Error in EnthalpyDifference::ObjectiveFunction : NaN detected\n";
				return ::MaxEnergy;
			}
		}
		if (x[0] < 0)
		{
			std::cerr << "Error in EnthalpyDifference::ObjectiveFunction : x[0]=" << x[0] << '\n';
			return ::MaxEnergy;
		}

		assert(grad==nullptr);
		AuxType * par =  reinterpret_cast<AuxType *> (aux);

		assert(par->pParam->NumParameters() == n);

		par->NumEvaluation++;
		std::memcpy(par->pParam->Parameters, x, sizeof(double)*par->pParam->NumParameters());
		IsotropicPairPotential * pot = par->pParam->GetPotential();

		double Pressure = GetPressure(par->pvpSolvers->at(0), par->TargetStructureVolume, *pot, *(par->pvpSolvers), par->pParam->MinDistance() );

		double result=EnthalpyDifference(*par->pParam, *par->pvpSolvers, par->TargetStructureVolume, Pressure)/(1.0+par->pParam->Penalty());

		if(par->AlsoConsiderElastic)
		{
			if(result>0)
			{
				logfile<<"Target not globally stable!\n";
				result+= ::MaxEnergy;
			}
			else
			{
				Configuration tar=par->pvpSolvers->at(0)->GetStructure();
				RelaxStructure(tar, *pot, Pressure, par->pParam->MinDistance());
				double LowestEigen=LowestEigenValueOfDeformations(tar, *pot, Pressure);
				if(LowestEigen<0)
				{
					//logfile<<"ERROR: Target not locally stable, return MaxEnergy!\n";
					//result= ::MaxEnergy;
					logfile<<"ERROR: Target not locally stable, return zero!\n";
					result= 0.0;
				}
				else
					result*=LowestEigen;
			}
		}

		logfile<<"at evaluation"<<par->NumEvaluation<<", Parameters are:";
		for(size_t i=0; i<n; i++)
			logfile<<x[i]<<" \t";
		logfile<<"\nEnthalpyDiff::ObjFunc return:"<<result<<", Pressure="<<Pressure<<"\n\n";

		delete pot;
		return result;
	}

	double Optimize(const std::vector<LatticeSumSolver *> & vpSolvers, ParameterSet & param, double & Pressure, double TargetVolume, int SpecificStage)
	{
		AuxType aux;
		aux.NumEvaluation=0;
		aux.TargetStructureVolume=TargetVolume;
		aux.pParam = & param;
		aux.pvpSolvers = & vpSolvers;
		aux.AlsoConsiderElastic=false;

		size_t dimension = param.NumParameters();
		double * variable = new double[dimension];
		std::memcpy(variable, param.Parameters, param.NumParameters()*sizeof(double));
		double * LowerBounds = new double[dimension];
		double * UpperBounds = new double[dimension];
		double minf; // the minimum objective value, upon return

		for(int Stage=0; Stage<4; Stage++)
		{
			//stages:
			//0: Optimize so that Target is locally stable, use global optimization algorithm
			//1: Optimize so that Target is locally stable
			//2: Optimize so that Target is globally stable
			//3: Optimize to maximize dE*elasticEigen

			if(SpecificStage!=(-1) && SpecificStage!=Stage)
				continue;

			if(Stage==3)
			{
				Configuration tar=vpSolvers.at(0)->GetStructure();
				IsotropicPairPotential * pot = param.GetPotential();
				Pressure = GetPressure(vpSolvers.at(0), TargetVolume, *pot, vpSolvers, param.MinDistance());
				double LEigen=LowestEigenValueOfDeformations(tar, *pot, Pressure);
				delete pot;
				if(LEigen>0)
				{
					logfile<<"Another optimization, considering elastic\n";
					aux.AlsoConsiderElastic=true;
				}
				else
				{
					logfile<<"Lowest Eigenvalue of deformations="<<LEigen<<'\n';
					logfile<<"Target not stable yet, skip elastic optimization\n";
					delete [] variable;
					delete [] UpperBounds;
					delete [] LowerBounds;
					return ::MaxEnergy;
					break;
				}
			}
			//debug temp
			//if (Stage == 2)
			//{
			//	std::cout << "at stage 2, parameters are:\n";
			//	param.Write(std::cout);
			//	std::cout << '\n';
			//}

			nlopt_opt opt;
			if(Stage!=0)
				opt = nlopt_create(NLOPT_LN_COBYLA, dimension); 
			else
				opt = nlopt_create(NLOPT_GN_CRS2_LM, dimension);
			nlopt_result status;
			if(Stage>1)
			{
				status=nlopt_set_min_objective(opt, ObjectiveFunction, reinterpret_cast<void *>(&aux));
				assert(!(status<0));
				Configuration tar=vpSolvers.at(0)->GetStructure();
				DimensionType dim=tar.GetDimension();
				size_t nbr=tar.NumParticle();
				size_t NumInEqualityConstraints = dim*(nbr+dim)*2+(dim*(dim+1))/2+(nbr-1)*dim;
				std::vector<double> tols(NumInEqualityConstraints, 1e-3);
				status=nlopt_add_inequality_mconstraint(opt, NumInEqualityConstraints, InEquality_MultiConstraint, reinterpret_cast<void *>(&aux), &tols[0]);
				assert(!(status<0));
				status=nlopt_set_ftol_abs(opt, 1e-3);
				assert(!(status<0));
				status=nlopt_set_maxeval(opt, 20000);
				assert(!(status<0));
			}
			else
			{
				status=nlopt_set_min_objective(opt, EqualityConstraint, reinterpret_cast<void *>(&aux));
				assert(!(status<0));
				status=nlopt_set_stopval(opt, 1e-9);
				assert(!(status<0));
				status=nlopt_set_ftol_abs(opt, 1e-14);
				assert(!(status<0));
				status=nlopt_set_maxeval(opt, 2000);
				assert(!(status<0));
			}

			param.GetBound(LowerBounds, UpperBounds);
			status=nlopt_set_lower_bounds(opt, LowerBounds);
			assert(!(status<0));
			status=nlopt_set_upper_bounds(opt, UpperBounds);
			assert(!(status<0));
			try
			{
				status=nlopt_optimize(opt, variable, &minf);
			}
			catch (NotANumberFound a)
			{
				std::cerr<<"write parameterization to screen!\n";
				param.Write(std::cerr);
				std::cerr<<"Quitting\n";
				assert(false);
			}
			if(status<0 && status!=(-4))
				std::cerr<<"Nlopt error in EnthalpyDiffernece::Optimize, return value="<<status<<'\n';
			nlopt_destroy(opt);
			if(Stage==1 && minf>1e-3)
			{
				logfile<<"Local Stability not achieved, skip next steps!\n";
				delete [] variable;
				delete [] UpperBounds;
				delete [] LowerBounds;
				return ::MaxEnergy;
			}
			if(Stage==2 && minf>(-1e-9))
			{
				logfile<<"Global Stability not achieved, skip next steps!\n";
				delete [] variable;
				delete [] UpperBounds;
				delete [] LowerBounds;
				return ::MaxEnergy;
			}
		}

		//copy results back
		std::memcpy(param.Parameters, variable, param.NumParameters()*sizeof(double));

		if(SpecificStage==-1)
		{
			//double check
			//check if everything is correct after we have performed all optimization steps.
			double Repeat=ObjectiveFunction(dimension, variable, nullptr, reinterpret_cast<void *>(&aux));
			if( std::abs(Repeat-minf)/(std::abs(minf)+1) > 1e-5 )
			{
				std::cerr<<"Error in EnthalpyDifference::Optimize : Optimization result cannot be repeated\n";
				std::cerr<<minf<<" \t"<<Repeat<<'\n';
				::Plots(param, aux.pvpSolvers, "Error");
				minf= ::MaxEnergy;
			}
			Configuration tar=vpSolvers.at(0)->GetStructure();
			IsotropicPairPotential * pot = param.GetPotential();
			Pressure = GetPressure(vpSolvers.at(0), TargetVolume, *pot, vpSolvers, param.MinDistance());
			bool LocalStability=LocallyStable(tar, *pot, Pressure, param.MinDistance());
			if(LocalStability==false)
			{
				logfile<<"ERROR: Target not locally stable, discard this optimization!\n";
				//std::cerr<<"ERROR: Target not locally stable, discard this optimization!\n";
				minf= ::MaxEnergy;
			}
			DimensionType dim=tar.GetDimension();
			size_t nbr=tar.NumParticle();
			size_t NumInEqualityConstraints = dim*(nbr+dim)*2+(dim*(dim+1))/2+(nbr-1)*dim;
			std::vector<double> mconsresults(NumInEqualityConstraints, 0.0);
			InEquality_MultiConstraint(NumInEqualityConstraints, &mconsresults[0], dimension, variable, nullptr, reinterpret_cast<void *>(&aux));
			for(auto iter=mconsresults.begin(); iter!=mconsresults.end(); iter++)
			{
				if((*iter)>0.01)
				{
					//std::cerr<<"Error in EnthalpyDifference::Optimize : Optimization result violate constraints\n";
					logfile<<"Error in EnthalpyDifference::Optimize : Optimization result violate constraints\n";
					//std::cerr<<(*iter)<<'\n';
					//::Plots(param, aux.pvpSolvers, "Error");
					minf= ::MaxEnergy;
				}
			}
			delete pot;

		}
		delete [] variable;
		delete [] UpperBounds;
		delete [] LowerBounds;

		return minf;
	}
};



//a class to record all annoying Competitors
class CompetitorRecord
{
private:
	struct record
	{
		size_t HitTime;
		FromStructure solver;
		record(Configuration & list) : HitTime(1), solver(list)
		{
		}
	};
	std::vector<record> records;
public:
	void RegisterStructure(Configuration & list)
	{
		this->records.push_back(record(list));
		for(auto iter=this->records.begin(); iter!=this->records.end()-1; iter++)
		{
			if(SameSolver(&iter->solver, &this->records.back().solver))
			{
				iter->HitTime++;
				this->records.pop_back();
				break;
			}
		}
	}
	void PublishRecords(std::ostream & out)
	{
		std::sort(this->records.begin(), this->records.end(), [](const record & left, const record & right)->bool { return left.HitTime>right.HitTime; });
		for(auto iter=this->records.begin(); iter!=this->records.end(); iter++)
		{
			out<<"===========================================================================\n";
			out<<"A Competitor which has been hit "<<iter->HitTime<<" times.\nTheta Series:\n";
			iter->solver.UpdateTerms(5);
			for(auto iter2=iter->solver.Terms.begin(); iter2!=iter->solver.Terms.end(); iter2++)
				out<<iter2->number<<" \t"<<iter2->distance<<'\n';
			out<<".pos file:\n";
			Configuration stru=iter->solver.GetStructure();
			::Output(out, stru);
			out<<'\n';
		}

		this->records.clear();
	}
};


struct ReOptimize_TempStruct
{
	ParameterSet * pparam;
	double result;
	double CenterPressure;
	double DeltaPressure;
};

bool ReOptimize(ParameterSet * & pInputParam, const std::vector<LatticeSumSolver *> & Solvers, size_t ReOptimizeTime, double TargetStructureVolume, double & result, size_t & NumOutput, double & CenterPressure, double & DeltaPressure, RandomGenerator & gen, std::string OutputPrefix)
{
#ifdef RecordCompetitors
	//add some code to record the lowest Competitor
	CompetitorRecord rec;
#endif
	std::vector<ReOptimize_TempStruct> successes;
	DeltaPressure=0.0;

optimize:
#pragma omp parallel default(shared) num_threads(OptimizationNumThreads)
	{
		ParameterSet * pMyParam=pInputParam->clone();
		double DeltaPressure=0.0, CenterPressure;
		std::vector<LatticeSumSolver *> mySolvers;
		for(auto iter=Solvers.begin(); iter!=Solvers.end(); iter++)
		{
			Configuration temp=(*iter)->GetStructure();
			mySolvers.push_back(new FromStructure(temp));
		}
#pragma omp for schedule(dynamic)
		for(signed long i=0; i<ReOptimizeTime; i++)
		{
			logfile<<"Optimization trial number:"<<i<<'\n';
			std::cout<<'-';
			if(i>1)
				pMyParam->Randomize(gen);
			result=EnthalpyDifference::Optimize(mySolvers, *pMyParam, CenterPressure, TargetStructureVolume, (i==0)?3:(-1));

			if(result<-1e-10)
			{
				ReOptimize_TempStruct temp;
				temp.CenterPressure=CenterPressure;
				temp.DeltaPressure=DeltaPressure;
				temp.pparam=pMyParam->clone();
				temp.result=result;
#pragma omp critical
				{
					successes.push_back(temp);
					logfile<<"A successful result:";
					//std::cout<<"A successful result:"<<result<<'\n';
					pMyParam->Write(logfile);
					logfile<<"Continue optimization!\n";
				}
			}
			else
			{
				logfile<<"A unsuccessful result:";
				pMyParam->Write(logfile);
				logfile<<"Continue optimization!\n";

#ifdef RecordCompetitors
				//add some code to record the lowest Competitor
				Configuration temp=Solvers[0]->GetStructure();
				EnthalpyDifference::EnthalpyDifference(*pInputParam, Solvers, TargetStructureVolume, CenterPressure, &temp);
				rec.RegisterStructure(temp);
#endif
			}
		}
		for(auto iter=mySolvers.begin(); iter!=mySolvers.end(); iter++)
			delete *iter;
		delete pMyParam;
	}
	if(successes.empty())
	{
		bool evolveresult=pInputParam->Evolve();
		if(evolveresult==false)
		{
			std::cout<<"Can't evolve parameterizaiton any more, stop optimization!\n";
			logfile<<"Can't evolve parameterizaiton any more, stop optimization!\n";
			return false;
		}
		else
		{
			std::cout<<"Parameterization Evolved, continue optimization!\n";
			logfile<<"Parameterization Evolved, continue optimization!\n";

#ifdef RecordCompetitors
			std::fstream CompetitorLog("CompetitorLog.txt", std::fstream::app|std::fstream::out);
			CompetitorLog<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
			CompetitorLog<<"For parameterization:";
			pInputParam->Write(CompetitorLog);
			CompetitorLog<<"Encountered Competitors are:\n";
			rec.PublishRecords(CompetitorLog);
#endif
			goto optimize;
		}
	}
	else
	{
		std::nth_element(successes.begin(), successes.begin(), successes.end(), [](const ReOptimize_TempStruct & left, const ReOptimize_TempStruct & right)->bool{return left.result<right.result;});
		delete pInputParam;
		pInputParam=successes.begin()->pparam;
		CenterPressure=successes.begin()->CenterPressure;
		DeltaPressure=successes.begin()->DeltaPressure;
		result=successes.begin()->result;
		for(auto iter=successes.begin()+1; iter!=successes.end(); iter++)
		{
			delete iter->pparam;
		}

		std::cout<<'\n';
		std::cout<<"Lowest ObjectiveFunction:"<<result<<"Parameters:\n";
		pInputParam->Write(std::cout);
		std::cout<<"Write Plot files as Optimization_"<<NumOutput<<"...";
		logfile<<"Lowest ObjectiveFunction:"<<result<<"Parameters:\n";
		pInputParam->Write(logfile);
		logfile<<"Write Plot files as Optimization_"<<NumOutput<<"...";
		std::cout.flush();
		std::stringstream stemp;
		stemp<<OutputPrefix<<"Optimization_"<<NumOutput;
		char PrefixBuffer[63];
		stemp.getline(PrefixBuffer, 62);
		::Plots(*pInputParam, &Solvers, PrefixBuffer);
		NumOutput++;
		std::cout<<"done\nOptimization Finish!";
		logfile<<"done\nOptimization Finish!";
		return true;
	}
}

/*
struct ReOptimize_TempStruct
{
ParameterSet * pparam;
double result;
double CenterPressure;
double DeltaPressure;
};
bool compare(ReOptimize_TempStruct left, ReOptimize_TempStruct right)
{
return left.result<right.result;
}

void ReOptimize(ParameterSet * & pInputParam, const std::vector<LatticeSumSolver *> & Solvers, size_t ReOptimizeTime, double TargetStructureVolume, double & result, size_t & NumOutput, double & CenterPressure, double & DeltaPressure, RandomGenerator & gen)
{
//pre-calculation of theta series
for(auto iter = Solvers.begin(); iter!= Solvers.end(); iter++)
{
double temp=1;
while( (*iter)->Terms.size()<2)
{
(*iter)->UpdateTerms(temp);
temp *= 1.5;
}
double newRc = pInputParam->MaxDistance()/pInputParam->MinDistance() * (*iter)->Terms[1].distance;
(*iter)->UpdateTerms(newRc);
}

std::vector<ReOptimize_TempStruct> instances;
instances.resize(ReOptimizeTime);

{
bool EvolveResult=true;
const ParameterSet * Parent = pInputParam;
while(EvolveResult)
{
for(size_t i=0; i<ReOptimizeTime; i++)
{
instances[i].pparam=nullptr;
instances[i].result= ::MaxEnergy;
}
bool successful=false;
for(size_t i=0; i<ReOptimizeTime; i++)
{
ParameterSet * temp = Parent->clone();
instances[i].pparam=temp;
if(i!=0)
instances[i].pparam->Randomize(gen);

nlopt_srand(static_cast<long>(gen.RandomDouble()*65535));

instances[i].result=EnthalpyDifference::Optimize(Solvers, *temp, instances[i].CenterPressure, TargetStructureVolume);
instances[i].DeltaPressure=0;

logfile<<"-----------------------------------------------------------\n";
logfile<<"a result: Lowest Pressure:"<<instances[i].CenterPressure-instances[i].DeltaPressure<<", Highest Pressure"<<instances[i].CenterPressure+instances[i].DeltaPressure<<'\n';
temp->Write(logfile);
if(instances[i].result<-0.1)
successful=true;
if(instances[i].result<-0.01 && i>ReOptimizeTime/15)
successful=true;
if(successful)
break;
}

std::sort(instances.begin(), instances.end(), compare);

std::cout<<"Lowest ObjectiveFunction:"<<instances[0].result<<"Parameters:\n";
instances[0].pparam->Write(std::cout);
logfile<<"Lowest ObjectiveFunction:"<<instances[0].result<<"Parameters:\n";
instances[0].pparam->Write(logfile);
std::stringstream stemp;
stemp<<"Optimization_"<<NumOutput;
char PrefixBuffer[32];
stemp.getline(PrefixBuffer, 31);
std::cout<<"Write Plot files as Optimization_"<<NumOutput<<"...";
logfile<<"Write Plot files as Optimization_"<<NumOutput<<"...";
std::cout.flush();
::Plots(*instances[0].pparam, &Solvers, PrefixBuffer);
NumOutput++;
std::cout<<"done\n";
logfile<<"done\n";
if(instances[0].result<0)
{
std::cout<<"Parameters found, stop evolution\n";
logfile<<"Parameters found, stop evolution\n";
delete pInputParam;
pInputParam = instances[0].pparam;
result=instances[0].result;
CenterPressure=instances[0].CenterPressure;
DeltaPressure=instances[0].DeltaPressure;
EvolveResult=false;
}
else
{
EvolveResult=instances[0].pparam->Evolve();
if(EvolveResult==false)
{
std::cout<<"Can't evolve the Parameterization any more\n";
logfile<<"Can't evolve the Parameterization any more\n";
result=instances[0].result;
delete pInputParam;
pInputParam = instances[0].pparam;
CenterPressure=instances[0].CenterPressure;
DeltaPressure=instances[0].DeltaPressure;
}
}
if(EvolveResult)
{
Parent = instances[0].pparam->clone();
std::cout<<"Parameterization evolved. Continue optimization\n";
logfile<<"Parameterization evolved. Continue optimization\n";
}
for(size_t i=0; i<ReOptimizeTime; i++)
if(instances[i].pparam != pInputParam && instances[i].pparam != nullptr)
{
delete instances[i].pparam;
instances[i].pparam=nullptr;
}
}
}
}
*/
