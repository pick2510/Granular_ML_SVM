#include "StructureOptimization.h"

#define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <nlopt.h>

const double MaxVolume=10000000;
const double MinVolume=0.01;
const double MaxLength2=4000000;
const double MinVolumeCoeff = 0.0;//structure is invalid if Volume < MinVolumeCoeff*|a1|*|a2|*...*|ad|
size_t StructureOptimization_SaveConfigurationInterval=std::numeric_limits<size_t>::max();
bool StructureOptimization_SaveNumberedConfiguration = false;

bool StructureOptimization_SaveEquiEnergyConfiguration = false;
double StructureOptimization_SaveEquiEnergyConfigurationEnergyInterval=std::numeric_limits<double>::max();
std::string StructureOptimization_SaveEquiEnergyConfigurationName="";

//codes related to elastic constants
//\epsilon_{ij}=C_{ijkl}*e_{kl}
void PrintAllElasticConstants(std::ostream & out, Configuration stru, Potential & pot, double Pressure, bool InfTime, double * presult)
{
	DimensionType dim = stru.GetDimension();

	//RelaxStructure(stru, pot, Pressure, 0.0);

	for (DimensionType i = 0; i<dim; i++)
		for (DimensionType j = 0; j<dim; j++)
			for (DimensionType k = 0; k<dim; k++)
				for (DimensionType l = 0; l < dim; l++)
				{
					double c = ElasticConstant(pot, stru, i, j, k, l, Pressure, InfTime);
					out << "C" << i << j << k << l << "=" << c << '\n';
					if (presult != nullptr)
						presult[i*dim*dim*dim + j*dim*dim + k*dim + l] = c;
				}
}

double EnergyDerivativeToDeformation(Potential & pot, Configuration structure, DimensionType i, DimensionType j, double epsilon2, double Pressure, GeometryVector * newbasis)
{
	DimensionType dim=structure.GetDimension(); 
	try
	{
		pot.SetConfiguration(structure);
		double EnergyDerivative = 0;
		std::vector<double> grad;
		grad.resize(dim*dim);
		pot.EnergyDerivativeToBasisVectors(&grad[0], Pressure);
		for (DimensionType t = 0; t < dim; t++)
		{
			GeometryVector temp = structure.GetBasisVector(t);
			EnergyDerivative += grad[t*dim + i] * temp.x[j];
		}
		return EnergyDerivative;
		//return grad[j*dim + i];
	}
	catch (EnergyDerivativeToBasisVectors_NotImplemented & a)
	{
		pot.SetConfiguration(structure);
		double prevEnergy = pot.Energy() + Pressure*structure.PeriodicVolume();
		//second structure deformation to calculate force
		for (DimensionType t = 0; t < dim; t++)
		{
			newbasis[t] = structure.GetBasisVector(t);
			if (i == j)
				newbasis[t].x[i] += epsilon2*newbasis[t].x[j];
			else
			{
				double temp = newbasis[t].x[i];
				newbasis[t].x[i] += 0.5*epsilon2*newbasis[t].x[j];
				newbasis[t].x[j] += 0.5*epsilon2*temp;
			}
		}
		structure.ChangeBasisVector(&newbasis[0]);
		pot.SetConfiguration(structure);
		double afterEnergy = pot.Energy() + Pressure*structure.PeriodicVolume();
		//calculate result
		return (afterEnergy - prevEnergy)/ epsilon2;
	}
}
double ElasticConstant(Potential & pot, Configuration structure, DimensionType i, DimensionType j, DimensionType k, DimensionType l, double Pressure, bool InfiniteTime)
{
	const double epsilon1=1e-6, epsilon2=1e-9;
	DimensionType dim=structure.GetDimension(); 
	assert(i<dim);
	assert(j<dim);
	assert(k<dim);
	assert(l<dim);
	double Volume=structure.PeriodicVolume();
	GeometryVector newbasis[::MaxDimension];

	//deforme structure so that strain_{kl}=epsilon1
	for(DimensionType t=0; t<dim; t++)
	{
		newbasis[t]=structure.GetBasisVector(t);
		double temp=newbasis[t].x[k];
		newbasis[t].x[k]+=epsilon1*newbasis[t].x[l];
	}

	double PrevDerivative = EnergyDerivativeToDeformation(pot, structure, i, j, epsilon2, Pressure, newbasis);

	structure.ChangeBasisVector(&newbasis[0]);
	if(InfiniteTime)
	{
		if(Verbosity>=4)
			std::cout<<"Relax structure to calculate infinite-time elastic constant.\n";
		RelaxStructure_NLOPT(structure, pot, Pressure, 0, 0.0);
		RelaxStructure_MINOP_withoutPertubation(structure, pot, Pressure, 0, 0.0);
	}

	double AfterDerivative = EnergyDerivativeToDeformation(pot, structure, i, j, epsilon2, Pressure, newbasis);

	return (AfterDerivative - PrevDerivative) / Volume / epsilon1;
}


//EulerAngles should have size d(d-1)/2
gsl_matrix * GetRotation(std::vector<double> EulerAngles, int d)
{
	gsl_matrix * result = gsl_matrix_alloc(d, d);
	gsl_matrix * temp = gsl_matrix_alloc(d, d);
	gsl_matrix * temp2 = gsl_matrix_calloc(d, d);
	gsl_matrix_set_identity(result);
	std::vector<double>::iterator iter=EulerAngles.begin();
	for(int i=0; i<d; i++)
	{
		for(int j=i+1; j<d; j++)
		{
			gsl_matrix_set_identity(temp);
			double c=std::cos(*iter), s=std::sin(*iter);
			gsl_matrix_set(temp, i, i, c);
			gsl_matrix_set(temp, j, j, c);
			gsl_matrix_set(temp, i, j, s);
			gsl_matrix_set(temp, j, i, (-1.0)*s);
			iter++;

			gsl_matrix_swap(result, temp2);
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp2, temp, 0.0, result);
		}
	}
	gsl_matrix_free(temp);
	gsl_matrix_free(temp2);
	return result;
}

#include <gsl/gsl_matrix.h>
//it is possible to optimize this class by using symmetry of the elastic tensor, not implemented yet
class AllElasticConstants : public BigVector
{
public:
	DimensionType d;
	//relax the structure before calling this function
	AllElasticConstants(const Configuration & c, Potential & p, double Pressure) : BigVector(c.GetDimension()*c.GetDimension()*c.GetDimension()*c.GetDimension())
	{	
		this->d=c.GetDimension();
		for(DimensionType i=0; i<d; i++)
			for(DimensionType j=0; j<d; j++)
				for(DimensionType k=0; k<d; k++)
					for(DimensionType l=0; l<d; l++)
						this->x[i*d*d*d+j*d*d+k*d+l]=ElasticConstant(p, c, i, j, k, l, Pressure);
	}
	AllElasticConstants(DimensionType d) : BigVector(d*d*d*d)
	{
		this->d=d;
	}
	AllElasticConstants(const BigVector & src) : BigVector(src)
	{
		this->d=std::floor(0.5+std::pow(src.dim, 0.25));
	}
	AllElasticConstants(const AllElasticConstants & src) : BigVector(src), d(src.d)
	{
	}
	AllElasticConstants Rotate(gsl_matrix * rot)
	{
		assert(rot->size1==d);
		assert(rot->size2==d);
		AllElasticConstants result(*this);
		for(DimensionType i=0; i<d; i++)
			for(DimensionType j=0; j<d; j++)
				for(DimensionType k=0; k<d; k++)
					for(DimensionType l=0; l<d; l++)
					{
						result.x[i*d*d*d+j*d*d+k*d+l]=0;
						for(DimensionType p=0; p<d; p++)
							for(DimensionType q=0; q<d; q++)
								for(DimensionType r=0; r<d; r++)
									for(DimensionType s=0; s<d; s++)
										result.x[i*d*d*d+j*d*d+k*d+l]+=gsl_matrix_get(rot, i, p)*gsl_matrix_get(rot, j, q)*gsl_matrix_get(rot, k, r)*gsl_matrix_get(rot, l, s)*this->x[p*d*d*d+q*d*d+r*d+s];
					}

		return result;
	}
};

struct ElasticOptimizeStruct
{
	int d;//dimension of the configuration
	int n;//number of parameters in the elastic
	std::function<BigVector(const std::vector<double> & param)> GetElasticFunc;
	const AllElasticConstants * pElastic;
	size_t EvaluateCount;
};
double Elastic_Objective(unsigned n, const double *x, double *grad, void *data)
{
	ElasticOptimizeStruct * pData = reinterpret_cast<ElasticOptimizeStruct *>(data);
	assert(n==pData->d*(pData->d-1)/2+pData->n);
	assert(grad==nullptr);
	std::vector<double> EulerAngles, param;
	for(int i=0; i<pData->n; i++)
		param.push_back(x[i]);
	for(int i=pData->n; i<n; i++)
		EulerAngles.push_back(x[i]);
	gsl_matrix * rot = GetRotation(EulerAngles, pData->d);
	AllElasticConstants el=pData->GetElasticFunc(param);
	AllElasticConstants el2=el.Rotate(rot);
	gsl_matrix_free(rot);
	double result=(el2-*pData->pElastic).Modulus2()/pData->pElastic->Modulus2();
	if(Verbosity>=5 || ((Verbosity>=4) && (pData->EvaluateCount%100000==0)))
		std::cout<<"f="<<result<<'\n';
	pData->EvaluateCount++;
	return result;
}
void ElasticOptimize(std::function<BigVector(const std::vector<double> & param)> getElasticFunc, int nparam, std::vector<double> & param, std::vector<double> ubound, std::vector<double> lbound, const Configuration & c, Potential & pot, double pressure)
{
	AllElasticConstants target(c, pot, pressure);
	ElasticOptimizeStruct aux;
	aux.d=target.d;
	aux.n=nparam;
	aux.pElastic=&target;
	aux.GetElasticFunc=getElasticFunc;
	aux.EvaluateCount=0;

	int OptimizerDim=aux.d*(aux.d-1)/2+nparam;
	nlopt_opt opt=nlopt_create(NLOPT_LN_SBPLX, OptimizerDim);
	param.resize(OptimizerDim, 0.0);
	lbound.resize(OptimizerDim, (-1.0)*pi);
	ubound.resize(OptimizerDim, pi);
	nlopt_set_lower_bounds(opt, &lbound[0]);
	nlopt_set_upper_bounds(opt, &ubound[0]);
	nlopt_set_ftol_rel(opt, 1e-10);
	nlopt_set_min_objective(opt, Elastic_Objective, reinterpret_cast<void *>(&aux));

	double result;
	nlopt_optimize(opt, &param[0], &result);
	if(Verbosity>3)
		std::cout<<"Elastic Optimization complete, resulting distance="<<result<<'\n';
	else if(result > 1e-5)
		std::cout<<"Elastic Optimization complete, resulting distance="<<result<<", distance large, maybe model doesn't fit!\n";

	param.resize(nparam);
	nlopt_destroy(opt);
}

////////////////////////////////////////////////////////////////////////////
//codes related to RelaxStructure
class InvalidStructure : public std::exception
{
public:
	InvalidStructure()
	{
	}
	virtual const char * what() const throw()
	{
		return "Optimizer want to try an Invalid Structure!\n";
	}
};


namespace
{
	struct Aux
	{
		Potential * pPot;
		size_t NumParticle;
		DimensionType Dim;
		double * pPressure, *pMinDisntace;
		size_t NumEvaluation;
		Configuration currentConfig;
		int Switch;//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
		double RelativeCoordinatesRescale;
		size_t NumParam;
		double PrevObjective;

		//members regarding enabling neighbor list
		size_t ShortMoveCount;
		std::vector<double> PrevX0;
		Aux(Potential * pPot, const Configuration & orig) : currentConfig(orig), pPot(pPot)
		{
			ShortMoveCount = 0;
			PrevX0 = std::vector<double>(orig.NumParticle(), ::MaxDistance);
			pPot->SetConfiguration(orig);
			PrevObjective = pPot->Energy();
		}
	};
	void GetConfig(Aux * pAux, const double * x)
	{
		if (x[0] != x[0])
			throw NotANumberFound();
		DimensionType & dim = pAux->Dim;
		assert(dim <= ::MaxDimension);
		GeometryVector base[::MaxDimension];
		Configuration & list(pAux->currentConfig);
		if (pAux->Switch != 0)
		{
			for (DimensionType i = 0; i<dim; i++)
			{
				base[i].SetDimension(dim);
				for (DimensionType j = 0; j<dim; j++)
					base[i].x[j] = x[i*dim + j];

				if (base[i].Modulus2()> ::MaxLength2)
					throw InvalidStructure();
			}
			if (::Volume(&base[0], dim)<std::pow(*pAux->pMinDisntace, dim))
				throw InvalidStructure();
			if (::Volume(&base[0], dim)> ::MaxVolume)
				throw InvalidStructure();
			if (::Volume(&base[0], dim)< ::MinVolume)
				throw InvalidStructure();
			double cubeVolume = 1.0;
			for (DimensionType i = 0; i < dim; i++)
				cubeVolume *= base[i].Modulus2();
			cubeVolume = std::sqrt(cubeVolume);
			double Skewedness = ::Volume(&base[0], dim) / cubeVolume;
			//a very low skewedness decreases the efficiency of cell list
			if (Skewedness<MinVolumeCoeff)
				throw InvalidStructure();
			//debug temp
			//std::cout << "skewedness=" << Skewedness;
			list.ChangeBasisVector(&base[0]);
		}
		if (pAux->Switch == 2)
		{
			for (size_t i = 0; i<pAux->NumParticle; i++)
			{
				GeometryVector ParticleRelative(dim);
				for (DimensionType j = 0; j<dim; j++)
					ParticleRelative.x[j] = x[dim*dim + dim*i + j] / pAux->RelativeCoordinatesRescale;
				list.MoveParticle(i, ParticleRelative);
			}
		}
		if (pAux->Switch == 0)
		{
			for (size_t i = 0; i<pAux->NumParticle; i++)
			{
				GeometryVector ParticleRelative(dim);
				for (DimensionType j = 0; j<dim; j++)
					ParticleRelative.x[j] = x[dim*i + j] / pAux->RelativeCoordinatesRescale;
				list.MoveParticle(i, ParticleRelative);
			}
		}
		if (std::time(nullptr) > ::TimeLimit)
		{
			std::fstream ofile("StructureOptimization_Progress_TimeLimit.configuration", std::fstream::out | std::fstream::binary);
			list.WriteBinary(ofile);
			ofile.close();
			exit(0);
		}
		if ((pAux->NumEvaluation + 1) % (StructureOptimization_SaveConfigurationInterval) == 0 && StructureOptimization_SaveNumberedConfiguration)
		{
			//std::stringstream ss;
			//ss << "StructureOptimization_Progress_" << pAux->NumEvaluation<<".configuration";
			//std::fstream ofile(ss.str(), std::fstream::out | std::fstream::binary);
			//list.WriteBinary(ofile);
			ConfigurationPack pk("StructureOptimization_Progress");
			pk.AddConfig(list);
		}
		//if( (pAux->NumEvaluation+1)%(2*StructureOptimization_SaveConfigurationInterval)==StructureOptimization_SaveConfigurationInterval )
		//{
		//	std::fstream ofile("StructureOptimization_Progress_a.configuration", std::fstream::out | std::fstream::binary);
		//	list.WriteBinary(ofile);
		//}
		//else if( (pAux->NumEvaluation+1)%(2*StructureOptimization_SaveConfigurationInterval)==0 )
		//{
		//	std::fstream ofile("StructureOptimization_Progress_b.configuration", std::fstream::out | std::fstream::binary);
		//	list.WriteBinary(ofile);
		//}

		return;
	}

	double Objective(unsigned n, const double *x, double *grad, void *data)
	{
		Aux * pAux = reinterpret_cast<Aux *>(data);
		DimensionType & dim = pAux->Dim;
		assert(n == pAux->NumParam);

		try
		{
			GetConfig(pAux, x);
			Configuration & list = pAux->currentConfig;

			//decide if the configuration should be put to cell list mode or neighbor list mode
			bool shouldUseNeighborList;
			PairPotential * pPairPot = dynamic_cast<PairPotential *>(pAux->pPot);
			if (pPairPot == nullptr)
				shouldUseNeighborList = false;
			else
			{
				//double dx0 = 0;
				//std::vector<double> x0(list.NumParticle(), 0.0);
				//for (int i = 0; i < list.NumParticle(); i++)
				//{
				//	x0[i] = list.GetCartesianCoordinates(i).x[0];
				//	dx0 = std::max(dx0, x0[i] - pAux->PrevX0[i]);
				//}
				//std::swap(x0, pAux->PrevX0);
				//if (std::abs(dx0) < 0.02*pPairPot->Rcut)
				//	pAux->ShortMoveCount++;
				//else
				//	pAux->ShortMoveCount = 0;
				//shouldUseNeighborList = (pAux->ShortMoveCount > 3);

				shouldUseNeighborList = true;
			}
			if (shouldUseNeighborList && !list.UsingNeighborList())
				list.SwitchToNeighborListMode();
			//else if (!shouldUseNeighborList && list.UsingNeighborList())
			//	list.SwitchToCellListMode();

			double Volume = list.PeriodicVolume();
			if (Volume< std::pow(::LengthPrecision, static_cast<double>(dim)))
				return ::MaxEnergy;

			if (*pAux->pMinDisntace> LengthPrecision)
			{
				for (size_t i = 0; i<list.NumParticle(); i++)
				{
					//Configuration::particle * pa=list.GetParticle(i);
					bool TooNearNeighbor = false;
					list.IterateThroughNeighbors(i, *pAux->pMinDisntace, [&TooNearNeighbor, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) -> void
					{
						if (SourceAtom != i || LatticeShift.Modulus2() != 0)
							TooNearNeighbor = true;
					}, &TooNearNeighbor);
					if (TooNearNeighbor)
						return ::MaxEnergy;
				}
			}

			double Pressure = (*pAux->pPressure);

			pAux->pPot->SetConfiguration(list);
			double result;

			if (grad == nullptr)
			{
				result = pAux->pPot->Energy() + Pressure*Volume;
				if (Verbosity >= 5 || (Verbosity == 4 && pAux->NumEvaluation % 1000 == 999))
				{
					std::cout.precision(14);
					std::cout << "at time " << std::time(nullptr) - ProgramStart << ", n=" << pAux->NumEvaluation << ", f=" << result << ", V=" << list.PeriodicVolume();
				}
			}
			else
			{
				for (size_t i = 0; i<n; i++)
					grad[i] = 0.0;

				if (pAux->Switch == 2)
				{
					std::vector<GeometryVector> forces;
					result = pAux->pPot->AllForceAndEnergy(forces) + Pressure*Volume;
					for (size_t i = 0; i < list.NumParticle(); i++)
					{
						//Configuration::particle * pA=list.GetParticle(i);
						GeometryVector & Force = forces[i];

						//derivatives of atom coordinates
						for (DimensionType j = 0; j<dim; j++)
						{
							grad[dim*dim + i*dim + j] = (-1.0)*Force.Dot(list.GetBasisVector(j)) / pAux->RelativeCoordinatesRescale;
						}
					}
				}
				else if (pAux->Switch == 0)
				{
					std::vector<GeometryVector> forces;
					result = pAux->pPot->AllForceAndEnergy(forces) + Pressure*Volume;
					for (size_t i = 0; i<list.NumParticle(); i++)
					{
						GeometryVector & Force = forces[i];

						for (DimensionType j = 0; j<dim; j++)
						{
							grad[i*dim + j] = (-1.0)*Force.Dot(list.GetBasisVector(j)) / pAux->RelativeCoordinatesRescale;
						}
					}
				}
				if (pAux->Switch != 0)
				{
					pAux->pPot->EnergyDerivativeToBasisVectors(grad, Pressure);
				}

				if (Verbosity >= 5 || (Verbosity == 4 && pAux->NumEvaluation % 1000 == 999))
				{
					if (Verbosity > 8)
					{
						std::cout << ", g=";
						auto prec = std::cout.precision();
						std::cout.precision(3);
						for (size_t i = 0; i < n; i++)
							std::cout << grad[i] * pAux->RelativeCoordinatesRescale << ',';
						std::cout.precision(prec);
						std::cout << ", ";
					}
					std::cout.precision(14);
					std::cout << "at time " << std::time(nullptr) - ProgramStart << ", n=" << pAux->NumEvaluation << ", f=" << result << ", V=" << list.PeriodicVolume();
					double normg = 0.0;
					for (size_t i = 0; i<n; i++)
						normg += grad[i] * grad[i];
					std::cout << ", |g|=" << normg*pAux->RelativeCoordinatesRescale*pAux->RelativeCoordinatesRescale;
				}

			}

			if (Verbosity >= 5 || (Verbosity == 4 && pAux->NumEvaluation % 1000 == 999))
			{
				std::cout << '\n';
			}

			pAux->NumEvaluation++;

			if (StructureOptimization_SaveEquiEnergyConfiguration && (pAux->PrevObjective - result) > StructureOptimization_SaveEquiEnergyConfigurationEnergyInterval)
			{
				ConfigurationPack pk(StructureOptimization_SaveEquiEnergyConfigurationName);
				pk.AddConfig(list);
				pAux->PrevObjective = result;
			}

			return result;
		}
		catch (InvalidStructure a)
		{
			if (grad != nullptr)
				for (unsigned i = 0; i<n; i++)
					grad[i] = 0.0;
			return ::MaxEnergy;
		}
	}

	void Gradient(unsigned n, const double *x, double *grad, void *data)
	{
		Aux * pAux=reinterpret_cast<Aux *>(data);
		DimensionType & dim=pAux->Dim;
		assert(n==pAux->NumParam);

		try
		{
			GetConfig(pAux, x);
			Configuration & list = pAux->currentConfig;
			double Volume = list.PeriodicVolume();
			if(Volume< std::pow(::LengthPrecision, static_cast<double>(dim)))
				return ;

			if(*pAux->pMinDisntace> LengthPrecision)
			{
				for(size_t i=0; i<list.NumParticle(); i++)
				{
					//Configuration::particle * pa=list.GetParticle(i);
					bool TooNearNeighbor=false;
					list.IterateThroughNeighbors(i, *pAux->pMinDisntace, [&TooNearNeighbor, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) -> void 
					{
						if(SourceAtom!=i || LatticeShift.Modulus2()!=0)
							TooNearNeighbor=true;
					}, &TooNearNeighbor);
					if(TooNearNeighbor)
						return ;
				}
			}

			double Pressure = (*pAux->pPressure);

			pAux->pPot->SetConfiguration(list);
			//double result=pAux->pPot->Energy()+ Pressure*Volume;

			if(Verbosity>=5 ||  (Verbosity==4&&pAux->NumEvaluation%1000==999) )
			{
				std::cout.precision(14);
				std::cout<<"at time "<<std::time(nullptr)-ProgramStart<<", n="<<pAux->NumEvaluation;
			}

			if(grad != nullptr)
			{
				for(size_t i=0; i<n; i++)
					grad[i]=0.0;

				if(pAux->Switch==2)
				{
					std::vector<GeometryVector> forces;
					pAux->pPot->AllForce(forces);
					for(size_t i=0; i<list.NumParticle(); i++)
					{
						//Configuration::particle * pA=list.GetParticle(i);
						GeometryVector & Force=forces[i];

						//derivatives of atom coordinates
						for(DimensionType j=0; j<dim; j++)
						{
							grad[dim*dim+i*dim+j] = (-1.0)*Force.Dot(list.GetBasisVector(j))/pAux->RelativeCoordinatesRescale;
						}
					}
				}
				else if(pAux->Switch==0)
				{
					std::vector<GeometryVector> forces;
					pAux->pPot->AllForce(forces);
					for(size_t i=0; i<list.NumParticle(); i++)
					{
						GeometryVector & Force=forces[i];

						for(DimensionType j=0; j<dim; j++)
						{
							grad[i*dim+j] = (-1.0)*Force.Dot(list.GetBasisVector(j))/pAux->RelativeCoordinatesRescale;
						}
					}
				}
				if(pAux->Switch!=0)
				{
					pAux->pPot->EnergyDerivativeToBasisVectors(grad, Pressure);
				}

				if(Verbosity>=5 ||  (Verbosity==4&&pAux->NumEvaluation%1000==999) )
				{
					if (Verbosity > 8)
					{
						std::cout << ", g=";
						auto prec = std::cout.precision();
						std::cout.precision(3);
						for (size_t i = 0; i < n; i++)
							std::cout << grad[i] * pAux->RelativeCoordinatesRescale << ',';
						std::cout.precision(prec);
						std::cout << ", ";
					}
					double normg = 0.0;
					for(size_t i=0; i<n; i++)
						normg+=grad[i]*grad[i];
					std::cout<<"|g|="<<normg*pAux->RelativeCoordinatesRescale*pAux->RelativeCoordinatesRescale;
				}

			}

			if(Verbosity>=5 ||  (Verbosity==4&&pAux->NumEvaluation%1000==999) )
			{
				std::cout<<'\n';
			}

			pAux->NumEvaluation++;

			return ;
		}
		catch(InvalidStructure a)
		{
			if(grad!=nullptr)
				for(unsigned i=0; i<n; i++)
					grad[i]=0.0;
			return ;
		}
	}

	//convert the configuration into parameters
	double * GetParam(Configuration & List, int Switch, double rescale, size_t & NumParam)
	{
		DimensionType dim=List.GetDimension();
		if(Switch==2)
			NumParam=dim*dim+dim*List.NumParticle();
		else if(Switch==1)
			NumParam=dim*dim;
		else if(Switch==0)
			NumParam=dim*List.NumParticle();
		else
		{
			std::cerr<<"Error in StructureOptimization : GetParam() : value of Switch not supported!\n";
			return nullptr;
		}
		double * pParams = new double[NumParam];
		//double * LBounds = new double[NumParam];
		//double * UBounds = new double[NumParam];
		//double * StepSizes = new double[NumParam];
		//double * Tolerences = new double[NumParam];
		//set basis vector parameters
		if(Switch!=0)
		{
			for(DimensionType i=0; i<dim; i++)
			{
				GeometryVector vBase = List.GetBasisVector(i);
				double Length = std::sqrt(vBase.Modulus2());
				for(DimensionType j=0; j<dim; j++)
				{
					pParams[i*dim+j]=vBase.x[j];
					//LBounds[i*dim+j]=-HUGE_VAL;
					//UBounds[i*dim+j]=HUGE_VAL;
					//StepSizes[i*dim+j]=0.001*Length;
					//Tolerences[i*dim+j]=1e-10*Length;
				}
			}
			if(Switch==2)
			{
				for(size_t i=0; i<List.NumParticle(); i++)
				{
					for(DimensionType j=0; j<dim; j++)
					{
						//Configuration::particle * pA=List.GetParticle(i);
						pParams[dim*dim+i*dim+j] = List.GetRelativeCoordinates(i).x[j]*rescale;
						//LBounds[dim*dim+i*dim+j] = -HUGE_VAL;
						//UBounds[dim*dim+i*dim+j] = HUGE_VAL;
						//StepSizes[dim*dim+i*dim+j] = 0.001;
						//Tolerences[dim*dim+i*dim+j]=1e-10;
					}
				}
			}
		}
		else
		{
			for(size_t i=0; i<List.NumParticle(); i++)
			{
				for(DimensionType j=0; j<dim; j++)
				{
					//Configuration::particle * pA=List.GetParticle(i);
					pParams[i*dim+j] = List.GetRelativeCoordinates(i).x[j]*rescale;
					//LBounds[i*dim+j] = -HUGE_VAL;
					//UBounds[i*dim+j] = HUGE_VAL;
					//StepSizes[i*dim+j] = 0.001;
					//Tolerences[i*dim+j]=1e-10;
				}
			}
		}
		return pParams;
	}
};


void RelaxStructure_FIRE(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double dtMax)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	//parameters according to 
	size_t Nmin = 5;
	double finc = 1.1, fdec = 0.5, aStart = 0.1, fa = 0.99;

	pot.SetConfiguration(List);
	DimensionType dim = List.GetDimension();
	size_t NumParam;

	//do optimization
	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle = List.NumParticle();
	aux.pPressure = &Pressure;
	aux.pMinDisntace = &MinDistance;
	aux.NumEvaluation = 0;
	aux.Switch = Switch;
	aux.RelativeCoordinatesRescale = 1.0;
	void * pAux = reinterpret_cast<void *>(&aux);

	double * pParams = GetParam(List, Switch, 1.0, NumParam);
	aux.NumParam = NumParam;
	if (pParams == nullptr)
		return;

	double * gradient = new double[NumParam];
	double * v = new double[NumParam];
	std::memset(v, 0, sizeof(double)*NumParam);
	std::memset(gradient, 0, sizeof(double)*NumParam);
	//double f1 = Objective(NumParam, pParams, gradient, pAux);

	double dt=1e-3, a=aStart;
	size_t NumStepsPPositive = 0;
	double normF = 1e10;
	double normV;
	double lastE = ::MaxEnergy;
	size_t NumSameEnergy = 0;
	for (;;)
	{
		for (size_t i = 0; i < NumParam; i++)
			v[i] -= 0.5*gradient[i] * dt;
		for (size_t i = 0; i < NumParam; i++)
			pParams[i] += v[i] * dt;
		double E = Objective(NumParam, pParams, gradient, pAux);
		for (size_t i = 0; i < NumParam; i++)
			v[i] -= 0.5*gradient[i] * dt;

		double P = 0;
		for (size_t i = 0; i < NumParam; i++)
			P -= gradient[i] * v[i];
		double sumF = 0.0, sumV = 0.0;
		for (int i = 0; i < NumParam; i++)
		{
			sumF += gradient[i] * gradient[i];
			sumV += v[i] * v[i];
		}
		normF = std::sqrt(sumF);
		normV = std::sqrt(sumV);

		for (size_t i = 0; i < NumParam; i++)
			v[i] = (1 - a)*v[i] - a*gradient[i] / normF*normV;

		if (P > 0)
			NumStepsPPositive++;
		else
			NumStepsPPositive = 0;

		if (NumStepsPPositive > Nmin)
		{
			dt = std::min(dt*finc, dtMax);
			a = a*fa;
		}
		if (P < 0 || E>lastE)
		{
			dt = dt*fdec;
			std::memset(v, 0, sizeof(double)*NumParam);
			a = aStart;
		}

		if (MaxStep >= 2)
			MaxStep -= 2;
		else
		{
			//debug temp
			//std::cout << "break because of MaxStep\n";

			break;
		}

		if (E>=lastE)
		{
			NumSameEnergy++;
			if (NumSameEnergy > 50)
			{
				//debug temp
				//std::cout << "break because of same energy\n";

				break;
			}
		}
		else
		{
			NumSameEnergy = 0;
			lastE = E;
		}
	}
	delete[] gradient;
	delete[] v;

	try
	{
		GetConfig(&aux, pParams);
		List = aux.currentConfig;
		List.SwitchToCellListMode();
	}
	catch (InvalidStructure a)
	{
		std::cerr << "Warning in RelaxStructure_FIRE : caught InvalidStructure error, stop structure optimization!\n";
	};
	delete[] pParams;
}

void RelaxStructure_FIRE_Emin(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double dtMax, double Emin)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	//parameters according to 
	size_t Nmin = 5;
	double finc = 1.1, fdec = 0.5, aStart = 0.1, fa = 0.99;

	pot.SetConfiguration(List);
	DimensionType dim = List.GetDimension();
	size_t NumParam;

	//do optimization
	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle = List.NumParticle();
	aux.pPressure = &Pressure;
	aux.pMinDisntace = &MinDistance;
	aux.NumEvaluation = 0;
	aux.Switch = Switch;
	aux.RelativeCoordinatesRescale = 1.0;
	void * pAux = reinterpret_cast<void *>(&aux);

	double * pParams = GetParam(List, Switch, 1.0, NumParam);
	aux.NumParam = NumParam;
	if (pParams == nullptr)
		return;

	double * gradient = new double[NumParam];
	double * v = new double[NumParam];
	std::memset(v, 0, sizeof(double)*NumParam);
	std::memset(gradient, 0, sizeof(double)*NumParam);
	//double f1 = Objective(NumParam, pParams, gradient, pAux);

	double dt = 1e-8, a = aStart;
	size_t NumStepsPPositive = 0;
	double normF = 1e10;
	double normV;
	double lastE = ::MaxEnergy;
	size_t NumSameEnergy = 0;
	for (;;)
	{
		for (size_t i = 0; i < NumParam; i++)
			v[i] -= 0.5*gradient[i] * dt;
		for (size_t i = 0; i < NumParam; i++)
			pParams[i] += v[i] * dt;
		double E = Objective(NumParam, pParams, gradient, pAux);
		for (size_t i = 0; i < NumParam; i++)
			v[i] -= 0.5*gradient[i] * dt;

		double P = 0;
		for (size_t i = 0; i < NumParam; i++)
			P -= gradient[i] * v[i];
		double sumF = 0.0, sumV = 0.0;
		for (int i = 0; i < NumParam; i++)
		{
			sumF += gradient[i] * gradient[i];
			sumV += v[i] * v[i];
		}
		normF = std::sqrt(sumF);
		normV = std::sqrt(sumV);

		for (size_t i = 0; i < NumParam; i++)
			v[i] = (1 - a)*v[i] - a*gradient[i] / normF*normV;

		if (P > 0)
			NumStepsPPositive++;
		else
			NumStepsPPositive = 0;

		if (NumStepsPPositive > Nmin)
		{
			dt = std::min(dt*finc, dtMax);
			a = a*fa;
		}
		if (P < 0 || E>lastE)
		{
			dt = dt*fdec;
			std::memset(v, 0, sizeof(double)*NumParam);
			a = aStart;
		}

		if (MaxStep >= 2)
			MaxStep -= 2;
		else
		{
			//debug temp
			//std::cout << "break because of MaxStep\n";

			break;
		}

		if (E < Emin)
		{
			//debug temp
			//std::cout << "break because of Emin\n";

			break;
		}

		if (E >= lastE)
		{
			NumSameEnergy++;
			if (NumSameEnergy > 50)
			{
				//debug temp
				//std::cout << "break because of same energy\n";

				break;
			}
		}
		else
		{
			NumSameEnergy = 0;
			lastE = E;
		}
	}
	delete[] gradient;
	delete[] v;

	try
	{
		GetConfig(&aux, pParams);
		List = aux.currentConfig;
		List.SwitchToCellListMode();
	}
	catch (InvalidStructure a)
	{
		std::cerr << "Warning in RelaxStructure_FIRE : caught InvalidStructure error, stop structure optimization!\n";
	};
	delete[] pParams;
}



double RelaxStructure_LocalGradientDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double maxError, double defaultStepSize, double maxTime)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	pot.SetConfiguration(List);
	DimensionType dim = List.GetDimension();
	size_t NumParam;

	//do optimization
	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle = List.NumParticle();
	aux.pPressure = &Pressure;
	aux.pMinDisntace = &MinDistance;
	aux.NumEvaluation = 0;
	aux.Switch = Switch;
	aux.RelativeCoordinatesRescale = 1.0;
	void * pAux = reinterpret_cast<void *>(&aux);

	double * pParams = GetParam(List, Switch, 1.0, NumParam);
	aux.NumParam = NumParam;
	if (pParams == nullptr)
		return 0.0;
	double * gradient = new double[NumParam];
	double * gradient2 = new double[NumParam];
	double * newx = new double[NumParam];
	double StepSize = defaultStepSize;

	double f1 = Objective(NumParam, pParams, gradient, pAux);
	double norm;
	//normalize the gradient
	double sum = 0.0;
	for (int i = 0; i<NumParam; i++)
		sum += gradient[i] * gradient[i];
	norm = std::sqrt(sum);
	for (int i = 0; i<NumParam; i++)
		gradient[i] /= norm;

	double time = 0.0;

	while (StepSize>1e-15)
	{
		double timeIncrement = StepSize / norm;

		for (int i = 0; i < NumParam; i++)
			newx[i] = pParams[i] - gradient[i] * StepSize;
		double f2 = Objective(NumParam, newx, gradient2, pAux);


		double df = f1 - f2;
		double df2 = StepSize*norm;

		double error = std::abs(df - df2) / df2;
		//debug temp
		//std::cout << "error=" << error << std::endl;
		if (error > maxError)
		{
			StepSize *= 0.5;

			//debug temp
			//std::cout<<"Size decreased to"<<StepSize<<'\n';
		}
		else
		{
			time += timeIncrement;
			std::memcpy(pParams, newx, sizeof(double)*NumParam);
			std::memcpy(gradient, gradient2, sizeof(double)*NumParam);
			//normalize the gradient
			double sum = 0.0;
			for (int i = 0; i<NumParam; i++)
				sum += gradient[i] * gradient[i];

			norm = std::sqrt(sum);
			for (int i = 0; i<NumParam; i++)
				gradient[i] /= norm;

			f1 = f2;
			if (error < 0.5*maxError)
			{
				StepSize *= 1.1;
				//debug temp
				//std::cout<<"Size increased to"<<StepSize<<'\n';
			}
		}

		if (MaxStep >= 2)
			MaxStep -= 2;
		else
			break;

		if (time > maxTime)
			break;

		//if(MaxStep%1000==0)
		//std::cout << "physical time spent in local gradient descent so far: " << time << std::endl;
	}
	delete[] gradient;
	delete[] gradient2;
	delete[] newx;

	try
	{
		GetConfig(&aux, pParams);
		List = aux.currentConfig;
		List.SwitchToCellListMode();
	}
	catch (InvalidStructure a)
	{
		std::cerr << "Warning in RelaxStructure_LocalGradientDescent : caught InvalidStructure error, stop structure optimization!\n";
	};
	delete [] pParams;

	//std::cout << "physical time spent in local gradient descent: " << time << std::endl;
	return time;
}

#include <memory>
void RelaxStructure_LocalGradientDescent_StepSnapshots(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double maxError, double stepInterval, std::string packName)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	std::fstream ofile(packName + "_timeLog.txt", std::fstream::out);

	std::unique_ptr<ConfigurationPack> ppk = nullptr;
	size_t nConfig = 0;

	pot.SetConfiguration(List);
	DimensionType dim = List.GetDimension();
	size_t NumParam;

	//do optimization
	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle = List.NumParticle();
	aux.pPressure = &Pressure;
	aux.pMinDisntace = &MinDistance;
	aux.NumEvaluation = 0;
	aux.Switch = Switch;
	aux.RelativeCoordinatesRescale = 1.0;
	void * pAux = reinterpret_cast<void *>(&aux);

	double * pParams = GetParam(List, Switch, 1.0, NumParam);
	aux.NumParam = NumParam;
	if (pParams == nullptr)
		return;

	double * gradient = new double[NumParam];
	double * gradient2 = new double[NumParam];
	double * newx = new double[NumParam];
	double StepSize = stepInterval;
	double totalStepLength = 0.0;

	double f1 = Objective(NumParam, pParams, gradient, pAux);
	double norm;
	//normalize the gradient
	double sum = 0.0;
	for (int i = 0; i<NumParam; i++)
		sum += gradient[i] * gradient[i];
	norm = std::sqrt(sum);
	for (int i = 0; i<NumParam; i++)
		gradient[i] /= norm;

	double time = 0.0;

	while (StepSize>1e-15)
	{
		double timeIncrement = StepSize / norm;

		for (int i = 0; i<NumParam; i++)
			newx[i] = pParams[i] - gradient[i] * StepSize;

		double f2 = Objective(NumParam, newx, gradient2, pAux);


		double df = f1 - f2;
		double df2 = StepSize*norm;

		double error = std::abs(df - df2) / df2;
		if (error > maxError)
		{
			StepSize *= 0.5;
		}
		else
		{
			//accept this move
			time += timeIncrement;
			totalStepLength += StepSize;
			if (std::floor(totalStepLength / stepInterval) > nConfig)
			{
				if (ppk == nullptr)
				{
					ppk.reset(new ConfigurationPack(packName));
				}
				ppk->AddConfig(aux.currentConfig);
				nConfig = ppk->NumConfig();

				ofile << time << std::endl;
				time = 0.0;
			}

			std::memcpy(pParams, newx, sizeof(double)*NumParam);
			std::memcpy(gradient, gradient2, sizeof(double)*NumParam);
			//normalize the gradient
			double sum = 0.0;
			for (int i = 0; i<NumParam; i++)
				sum += gradient[i] * gradient[i];

			norm = std::sqrt(sum);
			for (int i = 0; i<NumParam; i++)
				gradient[i] /= norm;

			f1 = f2;
			if (error < 0.5*maxError)
			{
				StepSize *= 1.1;
				if (StepSize > stepInterval)
					StepSize = stepInterval;
			}
		}

		if (MaxStep >= 2)
			MaxStep -= 2;
		else
			break;
	}
	delete[] gradient;
	delete[] gradient2;
	delete[] newx;

	try
	{
		GetConfig(&aux, pParams);
		List = aux.currentConfig;
		List.SwitchToCellListMode();
	}
	catch (InvalidStructure a)
	{
		std::cerr << "Warning in RelaxStructure_LocalGradientDescent : caught InvalidStructure error, stop structure optimization!\n";
	};
	delete[] pParams;
}



void RelaxStructure_NLOPT_inner(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, double rescale, size_t MaxStep, double fmin)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	pot.SetConfiguration(List);
	double E = pot.Energy();
	if(E<fmin)
		return;
	DimensionType dim=List.GetDimension();
	size_t NumParam;
	nlopt_result status;

	//do optimization
	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle=List.NumParticle();
	aux.pPressure= & Pressure;
	aux.pMinDisntace= & MinDistance;
	aux.NumEvaluation=0;
	aux.Switch=Switch;
	aux.RelativeCoordinatesRescale=rescale;
	nlopt_opt opt;

	double * pParams=GetParam(List, Switch, rescale, NumParam);
	aux.NumParam=NumParam;
	if(pParams==nullptr)
		return;

	opt = nlopt_create(NLOPT_LD_LBFGS, NumParam);
	nlopt_set_vector_storage(opt, 100);
	//nlopt_set_lower_bounds(opt, LBounds);
	//nlopt_set_upper_bounds(opt, UBounds);
	nlopt_set_min_objective(opt, Objective, reinterpret_cast<void *>(&aux));
	//nlopt_set_xtol_abs(opt, Tolerences);
	//nlopt_set_ftol_rel(opt, ftol_rel);
	//if(xtol_rel > 0)
	//	nlopt_set_xtol_rel(opt, xtol_rel);
	//nlopt_set_maxtime(opt, 300.0);
	//double TypicalLength = std::pow(List.PeriodicVolume()/List.NumParticle(), 1.0/List.GetDimension());
	//nlopt_set_initial_step1(opt, 0.1*TypicalLength);
	nlopt_set_stopval(opt, fmin);
	nlopt_set_maxeval(opt, MaxStep);
	nlopt_set_initial_step1(opt, 1e-4);

	double minf; 

	//std::cout<<"NLopt fvalue before:"<<Objective(NumParam, pParams, nullptr, reinterpret_cast<void *>(&aux))<<'\n';

	status=nlopt_optimize(opt, pParams, &minf);

	//debug temp
	//std::cout<<"NLopt return:"<<status<<'\n';
	//std::cout<<"NLopt fvalue after:"<<minf<<'\n';

	try
	{
		GetConfig(&aux, pParams);
		List = aux.currentConfig;
	}
	catch(InvalidStructure a)
	{
		//std::cerr<<"Warning in RelaxStructure_NLOPT_inner : caught InvalidStructure error, stop structure optimization!\n";
	};
	nlopt_destroy(opt);
	delete [] pParams;
}

void RelaxStructure_NLOPT_NoDerivative(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, double rescale, size_t MaxStep, double fmin)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	pot.SetConfiguration(List);
	if(pot.Energy()<fmin)
		return;
	DimensionType dim=List.GetDimension();
	size_t NumParam;
	nlopt_result status;

	//do optimization
	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle=List.NumParticle();
	aux.pPressure= & Pressure;
	aux.pMinDisntace= & MinDistance;
	aux.NumEvaluation=0;
	aux.Switch=Switch;
	aux.RelativeCoordinatesRescale=rescale;
	nlopt_opt opt;

	double * pParams=GetParam(List, Switch, rescale, NumParam);
	aux.NumParam=NumParam;
	if(pParams==nullptr)
		return;

	opt = nlopt_create(NLOPT_LN_SBPLX, NumParam);
	//nlopt_set_lower_bounds(opt, LBounds);
	//nlopt_set_upper_bounds(opt, UBounds);
	//nlopt_set_initial_step(opt, StepSizes);
	nlopt_set_min_objective(opt, Objective, reinterpret_cast<void *>(&aux));
	//nlopt_set_xtol_abs(opt, Tolerences);
	nlopt_set_ftol_rel(opt, 1e-15);
	nlopt_set_xtol_rel(opt, 1e-15);
	//nlopt_set_maxtime(opt, 300.0);
	nlopt_set_stopval(opt, fmin);
	nlopt_set_maxeval(opt, MaxStep);

	double minf; 

	//std::cout<<"NLopt fvalue before:"<<Objective(NumParam, pParams, nullptr, reinterpret_cast<void *>(&aux))<<'\n';

	status=nlopt_optimize(opt, pParams, &minf);

	//debug temp
	//std::cout<<"NLopt return:"<<status<<'\n';
	//std::cout<<"NLopt fvalue after:"<<minf<<'\n';

	try
	{
		GetConfig(&aux, pParams);
		List = aux.currentConfig;
		List.SwitchToCellListMode();
	}
	catch(InvalidStructure a)
	{
		//std::cerr<<"Warning in RelaxStructure_NLOPT_inner : caught InvalidStructure error, stop structure optimization!\n";
	};
	nlopt_destroy(opt);
	delete [] pParams;
	//delete [] LBounds;
	//delete [] UBounds;
	//delete [] StepSizes;
	//delete [] Tolerences;
}


void RelaxStructure_NLOPT(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	if(Verbosity>=3)
	{
		pot.SetConfiguration(List);
		double E = pot.Energy();
		double V = List.PeriodicVolume();
		std::cout<<"Initial E="<<E<<", V="<<V<<", H="<<E+Pressure*V<<'\n';
	}
	//RelaxStructure_NLOPT_inner(List, pot, Pressure, Switch, MinDistance, 10000, MaxStep, (-1.0)*MaxEnergy);
	RelaxStructure_NLOPT_inner(List, pot, Pressure, Switch, MinDistance, 1, MaxStep, (-1.0)*MaxEnergy);
	RelaxStructure_NLOPT_inner(List, pot, Pressure, Switch, MinDistance, 0.01, MaxStep, (-1.0)*MaxEnergy);
	if(Verbosity>=3)
	{
		pot.SetConfiguration(List);
		double E = pot.Energy();
		double V = List.PeriodicVolume();
		std::cout<<"Final E="<<E<<", V="<<V<<", H="<<E+Pressure*V<<'\n';
	}
	List.SwitchToCellListMode();
}
void RelaxStructure_NLOPT_Emin(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double Emin)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	//RelaxStructure_NLOPT_inner(List, pot, Pressure, Switch, MinDistance, 10000, MaxStep, Emin);
	//pot.SetConfiguration(List);
	//if (pot.Energy() < Emin)
	//{
	//	List.SwitchToCellListMode();
	//	return;
	//}
	//RelaxStructure_NLOPT_inner(List, pot, Pressure, Switch, MinDistance, 1, MaxStep, Emin);
	//pot.SetConfiguration(List);
	//if(pot.Energy()<Emin)
	//{
	//	List.SwitchToCellListMode();
	//	return;
	//}
	RelaxStructure_NLOPT_inner(List, pot, Pressure, Switch, MinDistance, 0.01, MaxStep, Emin);
	List.SwitchToCellListMode();
}


/////////////////////////////////////
// interface for gsl
#include <gsl/gsl_multimin.h>
double my_f (const gsl_vector  * x, void * params)
{
	Aux * pParam = reinterpret_cast<Aux *>(params);
	return Objective(pParam->NumParam, x->data, nullptr, params);
}
void my_df (const gsl_vector * x, void * params,gsl_vector * df)
{
	Aux * pParam = reinterpret_cast<Aux *>(params);
	Gradient(pParam->NumParam, x->data, df->data, params);
}
void my_fdf (const gsl_vector * x, void * params,double *f,gsl_vector * df)
{
	Aux * pParam = reinterpret_cast<Aux *>(params);
	(*f)=Objective(pParam->NumParam, x->data, df->data, params);
}
void RelaxStructure_SteepestDescent(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	size_t iter = 0;
	volatile int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	DimensionType dim=List.GetDimension();
	size_t NumParam;

	//do optimization
	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle=List.NumParticle();
	aux.pPressure= & Pressure;
	aux.pMinDisntace= & MinDistance;
	aux.NumEvaluation=0;
	aux.Switch=Switch;
	aux.RelativeCoordinatesRescale=1.0;

	double * pParams=GetParam(List, Switch, 1.0, NumParam);
	aux.NumParam=NumParam;
	if(pParams==nullptr)
		return;

	my_func.n = NumParam;
	my_func.f = my_f;
	my_func.df = my_df;
	my_func.fdf = my_fdf;
	my_func.params = reinterpret_cast<void *>(&aux);

	x = gsl_vector_alloc (NumParam);
	std::memcpy(x->data, pParams, sizeof(*pParams)*NumParam);


	T = gsl_multimin_fdfminimizer_steepest_descent;
	s = gsl_multimin_fdfminimizer_alloc (T, NumParam);
	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-1);


	//tol should be 1e-4 rather than 1e-1 for bfgs2 algorithm in gsl
	//T = gsl_multimin_fdfminimizer_vector_bfgs2;
	//s = gsl_multimin_fdfminimizer_alloc (T, NumParam);
	//gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

	double PrevF= ::MaxEnergy;
	size_t SameFCount=0;
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);

		if (status)
			break;

		//status = gsl_multimin_test_gradient (s->gradient, 1e-13);

		if(s->f==PrevF)
		{
			SameFCount++;
			if(SameFCount>20)
				break;
		}
		else
		{
			PrevF=s->f;
			SameFCount=0;
		}
	}
	while (status == GSL_CONTINUE && iter < MaxStep);
	//std::cout << "status=" << status << std::endl;

	try
	{
		GetConfig(&aux, s->x->data);
		List = aux.currentConfig;
		List.SwitchToCellListMode();
	}
	catch(InvalidStructure a)
	{
		std::cerr<<"Warning in RelaxStructure_NLOPT_inner : caught InvalidStructure error, stop structure optimization!\n";
	};
	delete [] pParams;
	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);
}


void RelaxStructure_ConjugateGradient(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double minG, double minF)//Switch==1: move Basis Vectors, Switch==0: move atoms, Switch==2:move both
{
	size_t iter = 0;
	volatile int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	DimensionType dim=List.GetDimension();
	size_t NumParam;

	//do optimization
	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle=List.NumParticle();
	aux.pPressure= & Pressure;
	aux.pMinDisntace= & MinDistance;
	aux.NumEvaluation=0;
	aux.Switch=Switch;
	aux.RelativeCoordinatesRescale=1.0;

	double * pParams=GetParam(List, Switch, 1.0, NumParam);
	aux.NumParam=NumParam;
	if(pParams==nullptr)
		return;

	my_func.n = NumParam;
	my_func.f = my_f;
	my_func.df = my_df;
	my_func.fdf = my_fdf;
	my_func.params = reinterpret_cast<void *>(&aux);

	x = gsl_vector_alloc (NumParam);
	std::memcpy(x->data, pParams, sizeof(*pParams)*NumParam);


	T = gsl_multimin_fdfminimizer_conjugate_pr;
	s = gsl_multimin_fdfminimizer_alloc (T, NumParam);
	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.00001, 1e-1);


	//tol should be 1e-4 rather than 1e-1 for bfgs2 algorithm in gsl
	//T = gsl_multimin_fdfminimizer_vector_bfgs2;
	//s = gsl_multimin_fdfminimizer_alloc (T, NumParam);
	//gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

	double PrevF= ::MaxEnergy;
	size_t SameFCount=0;
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);

		if (status)
			{
				//std::cout << "Note in RelaxStructure_ConjugateGradient : optimizer break. \n";
				break;
			}

		status = gsl_multimin_test_gradient (s->gradient, minG);

		if (s->f <= minF)
		{
			std::cout << "Note in RelaxStructure_ConjugateGradient : hit minF, exiting. \n";
			break;
		}

		if(s->f==PrevF)
		{
			SameFCount++;
			if (SameFCount > 20)
			{
				std::cout << "Note in RelaxStructure_ConjugateGradient : same f hit 20 times, exiting. \n";
				break;
			}
		}
		else
		{
			PrevF=s->f;
			SameFCount=0;
		}
	}
	while (status == GSL_CONTINUE && iter < MaxStep);

	try
	{
		GetConfig(&aux, s->x->data);
		List = aux.currentConfig;
		List.SwitchToCellListMode();
	}
	catch(InvalidStructure a)
	{
		std::cerr<<"Warning in RelaxStructure_NLOPT_inner : caught InvalidStructure error, stop structure optimization!\n";
	};
	delete [] pParams;
	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);
}














///////////////////////////////////////////////////////
//  RelaxStructure_MINOP
class OptimFunc{
public:

	Aux * pAux;

	// Functions required by the MINOP algorithm:
	int getDataLength() const;

	double evalF(const gsl_vector* input) const;
	void evalG(const gsl_vector* input, gsl_vector* output) const;

	void normalize(gsl_vector* input) const;

	OptimFunc()
	{
	}
	// Destructor, only non-pure virtual method
	~OptimFunc();
};


// Functions required by the MINOP algorithm:
int OptimFunc::getDataLength() const
{
	size_t NumParam;
	if(this->pAux->Switch==2)
		NumParam=this->pAux->Dim*this->pAux->Dim+this->pAux->Dim*this->pAux->NumParticle;
	else if(this->pAux->Switch==1)
		NumParam=this->pAux->Dim*this->pAux->Dim;
	else if(this->pAux->Switch==0)
		NumParam=this->pAux->Dim*this->pAux->NumParticle;

	return NumParam;
}

double OptimFunc::evalF(const gsl_vector* input) const
{
	return Objective(this->getDataLength(), input->data, nullptr, reinterpret_cast<void *>(this->pAux));
}
void OptimFunc::evalG(const gsl_vector* input, gsl_vector* output) const
{
	Gradient(this->getDataLength(), input->data, output->data, reinterpret_cast<void *>(this->pAux));
}

void OptimFunc::normalize(gsl_vector* input) const
{
	//should NOT normalize because we don't know the cell size and whether basis vectors are included in input data

	if(this->pAux->Switch==2)
	{
		double * data=input->data;
		for(size_t i=0; i<this->pAux->NumParticle; i++)
		{
			for(DimensionType j=0; j<this->pAux->Dim; j++)
			{
				double & c=data[this->pAux->Dim*this->pAux->Dim+i*this->pAux->Dim+j];
				c-=std::floor(c);
			}
		}
	}
	else if (this->pAux->Switch==0)
	{
		double * data=input->data;
		for(size_t i=0; i<this->pAux->NumParticle; i++)
		{
			for(DimensionType j=0; j<this->pAux->Dim; j++)
			{
				double & c=data[i*this->pAux->Dim+j];
				c-=std::floor(c);
			}
		}
	}
}

OptimFunc::~OptimFunc() 
{ 
	/* Nothing to do here */ 
}


void RelaxStructure_MINOP_withoutPertubation(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep, double Emin)
{
	DimensionType dim=List.GetDimension();
	size_t NumParam;
	double * ptemp=::GetParam(List, Switch, 1.0, NumParam);
	gsl_vector * vParams = gsl_vector_calloc(NumParam);
	gsl_vector * vResults = gsl_vector_calloc(NumParam);
	double * pParams = vParams->data;
	std::memcpy(pParams, ptemp, NumParam*sizeof(double));
	delete [] ptemp;
	//add a small distortion?
	/*
	RandomGenerator gen;
	for(size_t i=0; i<NumParam; i++)
	pParams[i]+=(gen.RandomDouble()-0.5)*1e-7;
	*/

	Aux aux(&pot, List);
	aux.Dim = dim;
	aux.NumParticle=List.NumParticle();
	aux.pPressure= & Pressure;
	aux.pMinDisntace= & MinDistance;
	aux.NumEvaluation=0;
	aux.Switch=Switch;
	aux.RelativeCoordinatesRescale=1.0;
	aux.NumParam=NumParam;
	OptimFunc opt;
	opt.pAux=&aux;

	void runMINOP79(const OptimFunc& function, 
		const gsl_vector* startPos, gsl_vector* resultPos,
		int maxIter, double goal, double minG, int verbose);

	try
	{
		runMINOP79(opt, vParams, vResults, MaxStep, Emin, 1e-13, Verbosity-6);
	}
	catch (NotANumberFound except)
	{
		std::cerr<<"found NaN in MinopStructureOptimization, stopping!\n";
		std::cerr<<"Output structure to stdcerr\n";
		std::cerr.precision(17);
		::Output(std::cerr, List);
		throw;
	}
	catch (std::bad_alloc a)
	{
		std::cerr<<"Warning in MinopStructureOptimization : caught std::bad_alloc error, stop structure optimization!\n";
	}

	try
	{
		GetConfig(&aux, vResults->data);
		List = aux.currentConfig;
		List.SwitchToCellListMode();
	}
	catch(InvalidStructure & a)
	{
		std::cerr<<"Warning in MinopStructureOptimization : caught InvalidStructure error, stop structure optimization!\n";
	}

	gsl_vector_free(vParams);
	gsl_vector_free(vResults);
}

void RelaxStructure_MINOP(Configuration & List, Potential & pot, double Pressure, int Switch, double MinDistance, size_t MaxStep)
{
	RelaxStructure_MINOP_withoutPertubation(List, pot, Pressure, Switch, MinDistance, MaxStep);
	return;
}

void runMINOP79(const OptimFunc& function, 
		const gsl_vector* startPos, gsl_vector* resultPos,
		int maxIter, double goal, double minG, int verbose)
		/*
		* maxIter : the number of iterations at which to stop
		* goal : if the function goes under it, stop
		* verbose : if >= 2, display the optimization progress
		*/
{
	// Scalar variables
	const int n = function.getDataLength();

	double stepSize = 1e-3;          // Called \Delta in the paper
	double minStepSize = 1e-15;
	double maxStepSize = 1;

	double f_x;            // func(v_x)
	double f_xa;           // func(v_xa)
	double norm_g;         // Needed quite often, so...


	// Vector variables (all of them should start with v_)
	gsl_vector* v_x = gsl_vector_calloc(n);     // position
	gsl_vector* v_xa = gsl_vector_calloc(n);    // new position
	gsl_vector* v_dx = gsl_vector_calloc(n);    // position difference
	gsl_vector* v_g = gsl_vector_calloc(n);     // gradient
	gsl_vector* v_ga = gsl_vector_calloc(n);    // new gradient
	gsl_vector* v_dg = gsl_vector_calloc(n);    // gradient difference
	gsl_vector* v_n = gsl_vector_calloc(n);     // -H*g
	gsl_vector* v_p = gsl_vector_calloc(n);
	gsl_vector* v_w = gsl_vector_calloc(n);     // temporary vector
	gsl_vector* v_Gdx = gsl_vector_calloc(n);
	gsl_vector* v_Hdg = gsl_vector_calloc(n);


	// Matrix variables (all of which should start with m_)
	// Both of these are symmetric, so only the lower triangular
	// part of them is actually used.
	gsl_matrix* m_G = gsl_matrix_calloc(n,n);   // approximate Hessian
	gsl_matrix* m_H = gsl_matrix_calloc(n,n);   // inverse approximate Hessian




	////// START
	gsl_vector_memcpy(v_x, startPos);

	f_x = function.evalF(v_x);
	function.evalG(v_x, v_g);
	norm_g = gsl_blas_dnrm2(v_g);


	gsl_matrix_set_zero(m_G);
	for(int i = 0; i < n; i++){
		gsl_matrix_set(m_G, i, i, norm_g / stepSize);
	}
	// H = G^-1
	gsl_matrix_set_zero(m_H);
	for(int i = 0; i < n; i++){
		gsl_matrix_set(m_H, i, i, stepSize / norm_g);
	}


	int iter;
	for(iter = 1; iter <= maxIter && f_x > goal && norm_g > minG &&
		stepSize >= minStepSize; iter++){

			// If verbose == 3, then print the information at every iteration.
			// If verbose == 2, then print it up to 100 times.
			//if(verbose >= 2 && (100*iter+50)/maxIter > (100*iter-50)/maxIter ||
			//	verbose >= 3){
			//		printf("[%6d]  step=%.2e  f(x)=%.8e  |g|=%.4e\n", iter, stepSize, f_x, norm_g);
			//}


			// Repeat until we actually get a decrease of the function value
			// The limit is there to avoid infinite loops
			while(stepSize >= 0.1 * minStepSize){

				////// STEP 1

				// n = -H*g
				gsl_blas_dsymv(CblasLower, -1.0, m_H, v_g, 0.0, v_n);

				double norm_n = gsl_blas_dnrm2(v_n);

				if(norm_n < stepSize){
					gsl_vector_memcpy(v_dx, v_n);
				} else { // if(norm_n < stepSize)
					double gGg = 0.0;
					double gHg = 0.0;

					// Since G and H are symmetric, the contribution from the
					// lower part is equal to the contribution to the upper part.
					for(int i = 0; i < n; i++){
						// Lower part, counts twice
						for(int j = 0; j < i; j++){
							gGg += 2.0 * gsl_vector_get(v_g, i) *
								gsl_matrix_get(m_G, i, j) *
								gsl_vector_get(v_g, j);
							gHg += 2.0 * gsl_vector_get(v_g, i) *
								gsl_matrix_get(m_H, i, j) *
								gsl_vector_get(v_g, j);
						}

						// Diagonal part, counts once
						gGg += gsl_vector_get(v_g, i) *
							gsl_matrix_get(m_G, i, i) *
							gsl_vector_get(v_g, i);
						gHg += gsl_vector_get(v_g, i) *
							gsl_matrix_get(m_H, i, i) *
							gsl_vector_get(v_g, i);
					}

					double c = pow(norm_g,4) / (gGg * gHg);
					double t = 0.2 + 0.8 * c;

					//when iter==1, t*norm_n should be equal to stepSize
					//add ||iter==1 to go to this part even if we have numerical imprecision
					if (t * norm_n <= stepSize || iter==1) {
						for(int i = 0; i < n; i++){
							// dx = (stepSize / ||n||) * n
							gsl_vector_set(v_dx, i, 
								(stepSize / norm_n) * gsl_vector_get(v_n, i));
						}
						double temp = gsl_vector_get(v_dx, 0);
						if(temp!=temp)
							throw NotANumberFound();
					} else { // if (t * norm_n <= stepSize) 

						////// STEP 2

						// n = t * n
						gsl_blas_dscal(t, v_n);

						// p = -(||g||^2 / g*G*g) * g
						{
							double factor = -pow(norm_g,2) / gGg;
							for(int i = 0; i < n; i++){
								gsl_vector_set(v_p, i,
									factor * gsl_vector_get(v_g, i));
							}
						}

						double norm_p = gsl_blas_dnrm2(v_p);
						if (norm_p >= stepSize) {
							// dx = -(delta / ||g||) * g;
							double factor = -stepSize / norm_g;
							for(int i = 0; i < n; i++){
								gsl_vector_set(v_dx, i,
									factor * gsl_vector_get(v_g, i));
							}
							double temp = gsl_vector_get(v_dx, 0);
							if(temp!=temp)
								throw NotANumberFound();
						} else { // if (gsl_blas_dnrm2(v_p) >= stepSize)

							////// STEP 3

							// w = n - p
							for(int i = 0; i < n; i++){
								gsl_vector_set(v_w, i,
									gsl_vector_get(v_n, i) -
									gsl_vector_get(v_p, i));
							}

							// theta = (stepSize^2 - ||p||^2) / 
							//	  (p*w + sqrt((p*w)^2 + ||w||^2 * (stepSize^2 - ||p||^2)));
							double pw = 0.0;
							gsl_blas_ddot(v_p, v_w, &pw);
							double norm_w = gsl_blas_dnrm2(v_w);

							double theta = (stepSize*stepSize - norm_p*norm_p) /
								(pw + sqrt(pw*pw + norm_w*norm_w * 
								(stepSize*stepSize - norm_p*norm_p)));
							if(theta!=theta)
								theta=0.0;

							// dx = p + theta * w;
							for(int i = 0; i < n; i++){
								gsl_vector_set(v_dx, i,
									gsl_vector_get(v_p, i) +
									theta * gsl_vector_get(v_w, i));
							}
							double temp = gsl_vector_get(v_dx, 0);
							if(temp!=temp)
								throw NotANumberFound();
						} // if (gsl_blas_dnrm2(v_p) >= stepSize)

					} // if (t * norm_n <= stepSize) 

					double temp = gsl_vector_get(v_dx, 0);
					if(temp!=temp)
						throw NotANumberFound();

				} // if(norm_n < stepSize)

				////// STEP 4

				// xa = x + dx
				for(int i = 0; i < n; i++){
					gsl_vector_set(v_xa, i,
						gsl_vector_get(v_x, i) + gsl_vector_get(v_dx, i));
				}

				f_xa = function.evalF(v_xa);


				if (f_xa <= f_x) {
					// We reduced the value of the function, so get out of the loop
					break;
				} else { // if (f_xa < f_x)
					stepSize /= 2;
				} // if (f_xa < f_x)
			} // while(stepSize)


			////// STEP 5

			function.evalG(v_xa, v_ga);
			for(int i = 0; i < n; i++){
				gsl_vector_set(v_dg, i,
					gsl_vector_get(v_ga, i) - gsl_vector_get(v_g, i));
			}


			// These are needed both to decide the new value of stepSize
			// and for the update of G and H.
			gsl_blas_dsymv(CblasLower, 1.0, m_G, v_dx, 0.0, v_Gdx);
			gsl_blas_dsymv(CblasLower, 1.0, m_H, v_dg, 0.0, v_Hdg);

			double gdx;
			gsl_blas_ddot(v_g, v_dx, &gdx);
			double gadx;
			gsl_blas_ddot(v_ga, v_dx, &gadx);
			double dxdg;
			gsl_blas_ddot(v_dx, v_dg, &dxdg);
			double dxGdx;
			gsl_blas_ddot(v_dx, v_Gdx, &dxGdx);
			double dgHdg;
			gsl_blas_ddot(v_dg, v_Hdg, &dgHdg);

			// Note: fa - f is normally negative, same for right-hand side
			if(f_xa - f_x > 0.1 * (gdx + 0.5 * dxGdx)){
				stepSize = 0.5 * gsl_blas_dnrm2(v_dx);
			} else { // if(f_xa - f_x > 0.1 * (gdx + 0.5 * dxGdx))
				// Using v_w, as a temporary vector, as its value is no longer
				// needed past step 3.
				// w = dg - Gdx
				for(int i = 0; i < n; i++){
					gsl_vector_set(v_w, i, gsl_vector_get(v_dg, i)
						- gsl_vector_get(v_Gdx, i));
				} // for(int i)

				if(gdx >= 2.0*gadx ||
					gsl_blas_dnrm2(v_w) <= 0.5 * gsl_blas_dnrm2(v_dg)) {
						stepSize = 2.0 * gsl_blas_dnrm2(v_dx);
				} else {
					stepSize = gsl_blas_dnrm2(v_dx);
				}
			} // if(fa - f > 0.1 * (gdx + 0.5 * dxGdx))

			/*  I believe this was added as a test
			stepSize *= 2;
			if(stepSize > maxStepSize){
			stepSize = maxStepSize;
			}
			*/

			// After finding a new position which reduces the value of the
			// objective, we need to update the Hessian and inverse Hessian
			// approximations.

			////// STEP 6

			f_x = f_xa;
			gsl_vector_memcpy(v_x, v_xa);
			gsl_vector_memcpy(v_g, v_ga);
			norm_g = gsl_blas_dnrm2(v_g);

			// b0 = dg*dx
			double b0;
			gsl_blas_ddot(v_dg, v_dx, &b0);


			////// STEP X: BFGS matrices update
			// Instead of following the Dennis and Mei paper,
			// I followed Kaufman's code in using the BFGS update.

			if(b0 >= 1e-30){
				// Finally, we can update the G and H matrices...

				// G* = G + dg*dg'/(dx'*dg) - (G*dx)*(G*dx)' / (dx'*G*dx)
				// H* = G*^-1
				// H* = H + ((dx'*dg + dg'*H*dg) / (dx'*dg)^2) * dx*dx'
				//        - ((H*dg)*dx' + dx*(H*dg)') / (dx'*dg)
				if(dxdg!=0.0)
					gsl_blas_dsyr(CblasLower, 1.0 / dxdg, v_dg, m_G);
				if(dxGdx!=0.0)
					gsl_blas_dsyr(CblasLower, -1.0 / dxGdx, v_Gdx, m_G);

				if(dxdg!=0.0)
					gsl_blas_dsyr(CblasLower, (dxdg + dgHdg) / (dxdg*dxdg), v_dx, m_H);
				if(dxdg!=0.0)
					gsl_blas_dsyr2(CblasLower, -1.0 / dxdg, v_Hdg, v_dx, m_H);
			} else { // if(b0 >= 1e-30)
				// b0 is too small, don't update G and H
				if(verbose >= 3){
					std::cout << "Small b0 : " << b0 << "\n";
				}
			} // if(b0 >= 1e-30)
	} // for(iter)


	if(verbose >= 1){
		std::cout << "Iterated " << iter << " times.\n";
		std::cout << "Reason for ending: ";
		if(iter > maxIter){
			std::cout << "reached the maximum number of iterations.\n";
		}
		if(f_x <= goal){
			std::cout << "reached the desired value for the function.\n";
		}
		if(norm_g <= minG){
			std::cout << "gradient became too small.\n";
		}
		if(stepSize < minStepSize){
			std::cout << "step size became too small.\n";
		}
	}

	// We are done, so let's write the answer in its proper variable
	gsl_vector_memcpy(resultPos, v_x);
	function.normalize(resultPos);

	// Freeing vectors
	gsl_vector_free(v_x);
	gsl_vector_free(v_xa);
	gsl_vector_free(v_dx);
	gsl_vector_free(v_g);
	gsl_vector_free(v_ga);
	gsl_vector_free(v_dg);
	gsl_vector_free(v_n);
	gsl_vector_free(v_p);
	gsl_vector_free(v_w);
	gsl_vector_free(v_Gdx);
	gsl_vector_free(v_Hdg);

	// Freeing matrices
	gsl_matrix_free(m_G);
	gsl_matrix_free(m_H);

}


#ifdef USE_EIGEN
#include <Eigen/Dense>
//if pEigenVectors is not nullptr, then fill it with eigenvectors
void GetHessianEigens(Configuration c, Potential * pPot, std::vector<double> & results, std::vector< std::vector<double> > * pEigenVectors)
{
	DimensionType dim = c.GetDimension();
	size_t Num = c.NumParticle();
	size_t dimTensor = dim*Num;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hessian(dimTensor, dimTensor);


	//I am disabling this special path for ShiftedCCPotentials because it interferes with target-S(k) potentials, which include perfect-glass potentials
	//ShiftedCCPotential * pCCPot = dynamic_cast<ShiftedCCPotential *>(pPot);

	//if (pCCPot != nullptr)
	//{
	//	//for stealth potentials
	//	pCCPot->SetConfiguration(c);
	//	pCCPot->GetRho();
	//	double volume = c.PeriodicVolume();

	//	//#pragma omp parallel for schedule(dynamic) num_threads(3)
	//	for (signed long i = 0; i < Num; i++)
	//	{
	//		for (size_t j = 0; j <= i; j++)
	//		{
	//			auto iter2 = pCCPot->RhoReal.begin();
	//			auto iter3 = pCCPot->RhoImag.begin();

	//			for (auto iter = pCCPot->constraints.begin(); iter != pCCPot->constraints.end(); iter++, iter2++, iter3++)
	//			{
	//				GeometryVector ri = c.GetCartesianCoordinates(i);
	//				GeometryVector rij = ri - c.GetCartesianCoordinates(j);
	//				double coeff = std::cos(iter->k.Dot(rij));
	//				if (i == j)
	//					coeff -= ((*iter2)*std::cos(iter->k.Dot(ri)) + (*iter3)*std::sin(iter->k.Dot(ri)));
	//				coeff = 2 * iter->V*coeff/volume;

	//				for (DimensionType di = 0; di < dim; di++)
	//					for (DimensionType dj = 0; dj < dim; dj++)
	//						(Hessian(i*dim + di, j*dim + dj)) += coeff*iter->k.x[di] * iter->k.x[dj];

	//				if (i != j)
	//					for (DimensionType di = 0; di < dim; di++)
	//						for (DimensionType dj = 0; dj < dim; dj++)
	//							(Hessian(j*dim + dj, i*dim + di)) += coeff*iter->k.x[di] * iter->k.x[dj];
	//			}
	//		}
	//	}

	//	//debug temp
	//	//for (signed long i = 0; i < dimTensor; i++)
	//	//	for (size_t j = 0; j <= i/dim; j++)
	//	//		for (DimensionType d = 0; d < dim; d++)
	//	//		{
	//	//			std::cout << "Hessian(" << i << ", " << j*dim + d << ") = " << Hessian(i, j*dim + d) << '\n';
	//	//			std::cout << "Hessian(" << j*dim + d << ", " << i << ") = " << Hessian(j*dim + d, i) << '\n';
	//	//		}

	//}
	//else
	{
		//other potentials
		const double Displacement = std::pow(c.PeriodicVolume() / c.NumParticle(), (1.0) / c.GetDimension())*2e-5;
		for (signed long i = 0; i < dimTensor; i++)
		{
			size_t np = i / dim; //the index of the particle to move

			GeometryVector r = c.GetCartesianCoordinates(np);
			GeometryVector backup = c.GetRelativeCoordinates(np);
			r.x[i%dim] += Displacement;
			r = c.CartesianCoord2RelativeCoord(r);
			c.MoveParticle(np, r);
			pPot->SetConfiguration(c);
			for (size_t j = 0; j <= np; j++)
			{
				GeometryVector f;
				pPot->Force(f, j);
				f.MultiplyFrom(-1.0 / Displacement);
				for (DimensionType d = 0; d < dim; d++)
				{
					Hessian(i, j*dim + d) = f.x[d];
					Hessian(j*dim + d, i) = f.x[d];

					//debug temp
					//std::cout << "Hessian(" << i << ", " << j*dim + d << ") = " << f.x[d] << '\n';
					//std::cout << "Hessian(" << j*dim+d << ", " << i << ") = " << f.x[d] << '\n';
				}
			}
			c.MoveParticle(np, backup);
		}
		//debug temp
		//for (signed long i = 0; i < dimTensor; i++)
		//	for (size_t j = 0; j <= i / dim; j++)
		//		for (DimensionType d = 0; d < dim; d++)
		//		{
		//			std::cout << "Hessian(" << i << ", " << j*dim + d << ") = " << Hessian(i, j*dim + d) << '\n';
		//			std::cout << "Hessian(" << j*dim + d << ", " << i << ") = " << Hessian(j*dim + d, i) << '\n';
		//		}
	}

	Eigen::SelfAdjointEigenSolver<decltype(Hessian)> Solver(Hessian, (pEigenVectors == nullptr) ? Eigen::EigenvaluesOnly : Eigen::ComputeEigenvectors);
	if (Solver.info() != Eigen::Success)
	{
		std::cout << "Error in GetHessianEigens : Eigen failed\n";
		return;
	}

	results.clear();
	auto eigval = Solver.eigenvalues();
	if (pEigenVectors != nullptr)
	{
		pEigenVectors->clear();
		Hessian = Solver.eigenvectors();
	}

	for (size_t i = 0; i<dimTensor; i++)
	{
		double dtemp = eigval(i);
		results.push_back(dtemp);
		if (pEigenVectors != nullptr)
		{
			pEigenVectors->push_back(std::vector<double>());
			for (size_t j = 0; j<dimTensor; j++)
				pEigenVectors->back().push_back(Hessian(j, i));
		}
	}
}

#include <MatOp/SparseSymShiftSolve.h>
#include <SymEigsShiftSolver.h> 
void LowestNHessianEigens(const Configuration & c, PairPotential & pot, size_t n, std::vector<double> & result2, std::vector< std::vector<double> > * pEigenVectors)
{
	int dim = c.GetDimension();
	int nbr = c.NumParticle();

	typedef Eigen::Triplet<double> T;
	std::vector<T> vT;

	for (int i = 0; i<nbr; i++)
	{
		AtomInfo a1 = c.GetCharacteristics(i);
		const Configuration * pc = &c;
		c.IterateThroughNeighbors(i, pot.Rcut,
			[&dim, &vT, &i, &pot, &a1, &pc, &nbr](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t srcAtom) ->void
		{
			if (shift.Modulus2()< ::LengthPrecision* ::LengthPrecision)
				return;
			int j = srcAtom;
			AtomInfo a2 = pc->GetCharacteristics(j);
			for (int m = 0; m < dim; m++)
			{
				for (int n = 0; n < dim; n++)
				{
					vT.push_back(T(i*dim + m, i*dim + n, (-1.0)*pot.SecondDerivative(shift, m, n, a1, a2)));
					if (j < nbr)
						vT.push_back(T(i*dim + m, j*dim + n, pot.SecondDerivative(shift, m, n, a1, a2)));
				}
			}
		}
		);
	}

	Eigen::SparseMatrix<double> mat(nbr*dim, nbr*dim);
	mat.setFromTriplets(vT.begin(), vT.end());
	vT.clear();

	Spectra::SparseSymShiftSolve<double> op(mat);
	Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<double> > eigs(&op, n, 3*n, 0.0);
	eigs.init();
	int nconv = eigs.compute(1000000, 1e-30);
	// Retrieve results
	if (eigs.info() == Spectra::SUCCESSFUL)
	{
		auto evalues = eigs.eigenvalues();

		result2.clear();
		for (size_t i = 0; i < evalues.size(); i++)
			result2.push_back((-1)*evalues[i]);

		if (pEigenVectors != nullptr)
		{
			pEigenVectors->clear();
			auto evs = eigs.eigenvectors();
			for (size_t i = 0; i < n; i++)
			{
				pEigenVectors->push_back(std::vector<double>());
				for (size_t j = 0; j < nbr*dim; j++)
					pEigenVectors->back().push_back(evs.coeff(j, i));
			}
		}
	}
}


//#include <unsupported/Eigen/ArpackSupport>
//void LowestNHessianEigens(const Configuration & c, PairPotential & pot, size_t n, std::vector<double> & result2, std::vector< std::vector<double> > * pEigenVectors)
//{
//	int dim = c.GetDimension();
//	int nbr = c.NumParticle();
//
//	typedef Eigen::Triplet<double> T;
//	std::vector<T> vT;
//
//	for (int i = 0; i<nbr; i++)
//	{
//		AtomInfo a1 = c.GetCharacteristics(i);
//		const Configuration * pc = &c;
//		c.IterateThroughNeighbors(i, pot.Rcut,
//			[&dim, &vT, &i, &pot, &a1, &pc, &nbr](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t srcAtom) ->void
//		{
//			if (shift.Modulus2()< ::LengthPrecision* ::LengthPrecision)
//				return;
//			int j = srcAtom;
//			AtomInfo a2 = pc->GetCharacteristics(j);
//			for (int m = 0; m < dim; m++)
//			{
//				for (int n = 0; n < dim; n++)
//				{
//					vT.push_back(T(i*dim + m, i*dim + n, (-1.0)*pot.SecondDerivative(shift, m, n, a1, a2)));
//					if (j < nbr)
//						vT.push_back(T(i*dim + m, j*dim + n, pot.SecondDerivative(shift, m, n, a1, a2)));
//				}
//			}
//		}
//		);
//	}
//
//	Eigen::SparseMatrix<double> mat(nbr*dim, nbr*dim);
//	mat.setFromTriplets(vT.begin(), vT.end());
//	vT.clear();
//
//	//Spectra::SparseSymShiftSolve<double> op(mat);
//	//Spectra::SymEigsShiftSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<double> > eigs(&op, n, 3 * n, 0.0);
//	//eigs.init();
//	//int nconv = eigs.compute(1000000, 1e-15);
//
//	typedef Eigen::SparseMatrix<double> SparseMat;
//	typedef Eigen::SimplicialLDLT<SparseMat> SparseChol;
//	typedef Eigen::ArpackGeneralizedSelfAdjointEigenSolver <SparseMat, SparseChol> Arpack;
//	Arpack arpack;
//
//	arpack.compute(mat, n, "SM");
//
//	// Retrieve results
//	{
//		auto evalues = arpack.eigenvalues();
//
//		result2.clear();
//		for (size_t i = 0; i < evalues.size(); i++)
//			result2.push_back((-1)*evalues[i]);
//
//		//if (pEigenVectors != nullptr)
//		//{
//		//	pEigenVectors->clear();
//		//	auto evs = eigs.eigenvectors();
//		//	for (size_t i = 0; i < n; i++)
//		//	{
//		//		pEigenVectors->push_back(std::vector<double>());
//		//		for (size_t j = 0; j < nbr*dim; j++)
//		//			pEigenVectors->back().push_back(evs.coeff(j, i));
//		//	}
//		//}
//	}
//}


#else
void GetHessianEigens(Configuration pConfig, Potential * pPot, std::vector<double> & results, std::vector< std::vector<double> > * pEigenVectors)
{
	std::cerr << "Error in GetHessianEigens : Eigen not enabled!\n";
}

#endif