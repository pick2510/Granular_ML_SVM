#ifndef COLLECTIVECOORDINATEPOTENTIAL_INCLUDED
#define COLLECTIVECOORDINATEPOTENTIAL_INCLUDED

#include "Potential.h"

//#define HAVETWOFIXEDPARTICLES


// wanted to write a class for Ewald summation. However, it seems to be hard to find a precise documentation
// that is more general than those online. Online descriptions usually only deal with 1/r^2 potential
// When doing it, it might be good to refer to Uche's 2004 collective coordinate paper
// Or refer to Sal's duality paper.

// \phi = 1/(2V) \sum_{\mathbf k} u(\mathbf k) | \rho(\mathbf k) |^2 + this->Shift
// compared to definition in uche2004constraints, it is shifted and scaled
#include <omp.h>
class ShiftedCCPotential : public Potential
{
protected:
	const Configuration * pConfig;
	double Volume;
public:
	struct KPointConstraint
	{
		GeometryVector k;//absolute k, not relative to basis vectors
		double V;
		KPointConstraint(const GeometryVector & k, double V) : k(k), V(V)
		{
		}
	};
	//users or derived classes should fill this vector with k-constraints
	std::vector<KPointConstraint> constraints;

	//this value is added to the energy. It is initialized as 0, but user can change it.
	double Shift;

	//this class will manage these members, users can also read it.
	std::vector<double> RhoReal, RhoImag;
	bool RhoValid;
	virtual void GetRho(void)
	{
		assert(this->pConfig != nullptr);
		if (this->RhoValid)
		{
			return;
		}
		this->Volume = this->pConfig->PeriodicVolume();

		signed long end = this->constraints.size();
		this->RhoReal.resize(end);
		this->RhoImag.resize(end);

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for (signed long i = 0; i<end; i++)
		{
			double tr = 0, ti = 0;//real and imag parts of result
#ifdef HAVETWOFIXEDPARTICLES
			tr+=2;
#endif
			for (size_t j = 0; j<pConfig->NumParticle(); j++)
			{
				GeometryVector c = pConfig->GetCartesianCoordinates(j);
				tr += std::cos(this->constraints[i].k.Dot(c));
				ti += std::sin(this->constraints[i].k.Dot(c));
			}
			this->RhoReal[i] = tr;
			this->RhoImag[i] = ti;
		}

		this->RhoValid = true;
	}

	ShiftedCCPotential(DimensionType dimension):Potential(dimension)
	{
		this->pConfig = nullptr;
		this->RhoValid = false;
		this->ParallelNumThread = 1;
		this->Shift = 0.0;
	}
	virtual ~ShiftedCCPotential()
	{
	}
	ShiftedCCPotential * clone()
	{
		return new ShiftedCCPotential(*this);
	}
	virtual double Energy()//return the total potential energy of the given configuration
	{
		if (this->RhoValid == false)
			this->GetRho();

		double result = 0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end = this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread = omp_get_thread_num();
#pragma omp for schedule(static)
			for (signed long i = 0; i<end; i++)
				threadresult[numthread] += this->constraints[i].V*(this->RhoReal[i] * this->RhoReal[i] + this->RhoImag[i] * this->RhoImag[i]);
		}
		for (auto iter = threadresult.begin(); iter != threadresult.end(); iter++)
			result += *iter;

		return result / (Volume)+this->Shift;
	}
	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		if (this->RhoValid == false)
			this->GetRho();

		result = GeometryVector(this->Dimension);

		signed long end = this->constraints.size();

		for (signed long j = 0; j<end; j++)
		{
			GeometryVector k = this->constraints[j].k;
			GeometryVector r = this->pConfig->GetCartesianCoordinates(i);
			result.MinusFrom((2.0)*this->constraints[j].V*(this->RhoImag[j] * std::cos(k.Dot(r)) - this->RhoReal[j] * std::sin(k.Dot(r)))*k);
		}

		result.MultiplyFrom(1.0 / Volume);
	}
	virtual void AllForce(std::vector<GeometryVector> & results)
	{
		if (this->RhoValid == false)
			this->GetRho();
		long end = this->pConfig->NumParticle();
		results.resize(end);
#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		for (long i = 0; i<end; i++)
			this->Force(results[i], i);
	}
	//Energy derivative respect to basis vectors
	//grad[n*dim+m] is the energy derivative of mth (coordinate) component of nth basis vector
	//virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	//{
	//	std::cerr<<"Warning : CollectiveCoordinatePotential::EnergyDerivativeToBasisVectors is not implemented!\n";
	//}
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->Dimension == Position.GetDimension());
		this->pConfig = &Position;
		this->RhoValid = false;
	}
};

//user need to set mt (Mtilde(k) of each component), 
class MultiShiftedCCPotential: public Potential
{
protected:
public:
	struct KPointConstraint
	{
		GeometryVector k;//absolute k, not relative to basis vectors
		double V;
		KPointConstraint(const GeometryVector & k, double V) : k(k), V(V)
		{
		}
	};
	std::vector<KPointConstraint> constraints;
	std::map<AtomInfo, std::vector<double>> mt;

	const Configuration * pConfig;
	std::vector<double> RhoReal, RhoImag;
	bool RhoValid;
	double Volume;
	void GetRho(void)
	{
		assert(this->pConfig != nullptr);
		if (this->RhoValid)
		{
			return;
		}
		this->Volume = this->pConfig->PeriodicVolume();

		signed long end = this->constraints.size();
		this->RhoReal.resize(end);
		this->RhoImag.resize(end);

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for (signed long i = 0; i<end; i++)
		{
			double tr = 0, ti = 0;//real and imag parts of result
#ifdef HAVETWOFIXEDPARTICLES
			tr += 2;
#endif
			for (size_t j = 0; j<pConfig->NumParticle(); j++)
			{
				GeometryVector c = pConfig->GetCartesianCoordinates(j);
				auto MtIter = mt.find(pConfig->GetCharacteristics(j));
				double m = MtIter->second.at(i);
				tr += std::cos(this->constraints[i].k.Dot(c))*m;
				ti += std::sin(this->constraints[i].k.Dot(c))*m;
			}
			this->RhoReal[i] = tr;
			this->RhoImag[i] = ti;
		}

		this->RhoValid = true;
	}

	//this value is added to the energy. It is initialized as 0, but user can change it.
	double Shift;

	MultiShiftedCCPotential(DimensionType dimension) :Potential(dimension)
	{
		this->pConfig = nullptr;
		this->RhoValid = false;
		this->ParallelNumThread = 1;
		this->Shift = 0.0;
	}
	virtual ~MultiShiftedCCPotential()
	{
	}
	MultiShiftedCCPotential * clone()
	{
		return new MultiShiftedCCPotential(*this);
	}
	virtual double Energy()//return the total potential energy of the given configuration
	{
		if (this->RhoValid == false)
			this->GetRho();

		double result = 0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end = this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread = omp_get_thread_num();
#pragma omp for schedule(static)
			for (signed long i = 0; i<end; i++)
				threadresult[numthread] += this->constraints[i].V*(this->RhoReal[i] * this->RhoReal[i] + this->RhoImag[i] * this->RhoImag[i]);
		}
		for (auto iter = threadresult.begin(); iter != threadresult.end(); iter++)
			result += *iter;

		return result / (Volume)+this->Shift;
	}
	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		if (this->RhoValid == false)
			this->GetRho();

		result = GeometryVector(this->Dimension);

		signed long end = this->constraints.size();
		auto MtIter = mt.find(this->pConfig->GetCharacteristics(i));
		for (signed long j = 0; j<end; j++)
		{
			GeometryVector k = this->constraints[j].k;
			GeometryVector r = this->pConfig->GetCartesianCoordinates(i);
			double m = MtIter->second.at(j);
			result.MinusFrom((2.0*m*m)*this->constraints[j].V*(this->RhoImag[j] * std::cos(k.Dot(r)) - this->RhoReal[j] * std::sin(k.Dot(r)))*k);
		}

		result.MultiplyFrom(1.0 / Volume);
	}
	virtual void AllForce(std::vector<GeometryVector> & results)
	{
		if (this->RhoValid == false)
			this->GetRho();
		long end = this->pConfig->NumParticle();
		results.resize(end);
		unsigned long tempNumThread = 1;
		std::swap(tempNumThread, this->ParallelNumThread);
#pragma omp parallel for num_threads(tempNumThread) default(shared)
		for (long i = 0; i<end; i++)
			this->Force(results[i], i);
		std::swap(tempNumThread, this->ParallelNumThread);
	}
	//Energy derivative respect to basis vectors
	//grad[n*dim+m] is the energy derivative of mth (coordinate) component of nth basis vector
	//virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	//{
	//	std::cerr<<"Warning : CollectiveCoordinatePotential::EnergyDerivativeToBasisVectors is not implemented!\n";
	//}
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->Dimension == Position.GetDimension());
		this->pConfig = &Position;
		this->RhoValid = false;
	}
};

//same as ShiftedCCPotential, but allows user to specify a weight for each particle
//rho=\sum_j weight_j * exp (i*k*r_j)
class ShiftedCCPotential_withWeight : public ShiftedCCPotential
{
protected:
	std::vector<double> weights;
public:
	ShiftedCCPotential_withWeight(const ShiftedCCPotential & Pot, const std::vector<double> & weights) : ShiftedCCPotential(Pot), weights(weights)
	{
	}
	virtual void GetRho(void)
	{
		assert(this->pConfig != nullptr);
		assert(this->pConfig->NumParticle() == weights.size());
		if (this->RhoValid)
		{
			return;
		}
		this->Volume = this->pConfig->PeriodicVolume();

		signed long end = this->constraints.size();
		this->RhoReal.resize(end);
		this->RhoImag.resize(end);

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for (signed long i = 0; i<end; i++)
		{
			double tr = 0, ti = 0;//real and imag parts of result
#ifdef HAVETWOFIXEDPARTICLES
			tr += 2;
#endif
			for (size_t j = 0; j<pConfig->NumParticle(); j++)
			{
				GeometryVector c = pConfig->GetCartesianCoordinates(j);
				tr += std::cos(this->constraints[i].k.Dot(c))*weights[j];
				ti += std::sin(this->constraints[i].k.Dot(c))*weights[j];
			}
			this->RhoReal[i] = tr;
			this->RhoImag[i] = ti;
		}
		this->RhoValid = true;
	}

	virtual void Force(GeometryVector & result, size_t i)//calculate the force for the ith Particle
	{
		if (this->RhoValid == false)
			this->GetRho();

		result = GeometryVector(this->Dimension);

		signed long end = this->constraints.size();

		for (signed long j = 0; j<end; j++)
		{
			GeometryVector k = this->constraints[j].k;
			GeometryVector r = this->pConfig->GetCartesianCoordinates(i);
			result.MinusFrom((2.0)*this->constraints[j].V*(this->RhoImag[j] * std::cos(k.Dot(r)) - this->RhoReal[j] * std::sin(k.Dot(r)))*k);
		}

		result.MultiplyFrom(weights[i] / Volume);
	}
	virtual void AllForce(std::vector<GeometryVector> & results)
	{
		if (this->RhoValid == false)
			this->GetRho();
		long end = this->pConfig->NumParticle();
		results.resize(end);
		unsigned long tempNumThread = 1;
		std::swap(tempNumThread, this->ParallelNumThread);
#pragma omp parallel for num_threads(tempNumThread) default(shared)
		for (long i = 0; i<end; i++)
			this->Force(results[i], i);
		std::swap(tempNumThread, this->ParallelNumThread);
	}
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->Dimension == Position.GetDimension());
		this->pConfig = &Position;
		this->RhoValid = false;
	}
};



// \phi = log ( 1/(2V) \sum_{\mathbf k} u(\mathbf k) | \rho(\mathbf k) |^2 )
class LogCCPotential : public ShiftedCCPotential
{
public:
	LogCCPotential( const ShiftedCCPotential & pot ) : ShiftedCCPotential(pot)
	{
	}
	virtual double Energy()
	{
		return std::log( ShiftedCCPotential::Energy() );
	}
	virtual void Force(GeometryVector & result, size_t i)
	{
		ShiftedCCPotential::Force(result, i);
		result.MultiplyFrom( 1.0/ShiftedCCPotential::Energy() );
	}
};

// \phi = \sum_{\mathbf k} V(\mathbf k) ( | \rho(\mathbf k) |^2/N - S0(\mathbf k) )^2 + this->Shift
class ShiftedCCPotential_S2 : public ShiftedCCPotential
{
public:
	std::vector<double> S0;
	size_t N;
	ShiftedCCPotential_S2(DimensionType d) : ShiftedCCPotential(d)
	{
	}
	virtual double Energy() //return the total potential energy of the given configuration
	{
		if(this->RhoValid==false)
			this->GetRho();

		double result=0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end=this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long i=0; i<end; i++)
			{
				double s=(this->RhoReal[i]*this->RhoReal[i] + this->RhoImag[i]*this->RhoImag[i])/N;
				double deltaS=(s-this->S0[i]);
				threadresult[numthread]+=this->constraints[i].V*deltaS*deltaS;
			}
		}
		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result+=*iter;

		result+=this->Shift;

		return result;
	}
	virtual void Force(GeometryVector & result, size_t i) //calculate the force for the ith Particle
	{
		if(this->RhoValid==false)
			this->GetRho();

		result=GeometryVector(this->Dimension);
		std::vector<GeometryVector> threadresult(this->ParallelNumThread, result);

		signed long end=this->constraints.size();

#pragma omp parallel num_threads(1) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long j=0; j<end; j++)
			{
				GeometryVector k=this->constraints[j].k;
				GeometryVector r=this->pConfig->GetCartesianCoordinates(i);
				double s=(this->RhoReal[j]*this->RhoReal[j] + this->RhoImag[j]*this->RhoImag[j])/N;
				double deltaS=(s-this->S0[j]);
				threadresult[numthread].MinusFrom((4.0/this->N)*this->constraints[j].V*deltaS*(this->RhoImag[j]*std::cos(k.Dot(r))-this->RhoReal[j]*std::sin(k.Dot(r)))*k);
			}
		}

		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result.AddFrom(*iter);
	}
	virtual void SetConfiguration(const Configuration & Position)
	{
		this->N=Position.NumParticle();
		this->ShiftedCCPotential::SetConfiguration(Position);
	}
};

// \phi = \sum_{\mathbf k} V(\mathbf k) ( <S(\mathbf k)> - S0(\mathbf k) )^2 + this->Shift
// where <...> is average over sub-configurations.
//treat the configuration as NSubConfig sub-configurations.
//particles 1 to NumParticle/NSubConfig is in the first sub-configuration
//<S(\mathbf k)> means S(\mathbf k) avaraged over all sub-configurations 
class ShiftedCCPotential_S2_MultiConfig : public ShiftedCCPotential
{
public:
	std::vector<double> S0;
	std::vector<double> Snow;//updated in GetRho
	size_t NSubConfig;//number of separate sub-configurations
	size_t Np;//number of particles per sub-configuration
	ShiftedCCPotential_S2_MultiConfig(DimensionType d, size_t NSubConfig) : ShiftedCCPotential(d), NSubConfig(NSubConfig)
	{
	}
	void GetRho_MultiConfig(void)
	{
		assert(this->pConfig!=nullptr);
		if(this->RhoValid)
		{
			return;
		}

		signed long end=this->constraints.size();
		this->RhoReal.resize(end*NSubConfig);
		this->RhoImag.resize(end*NSubConfig);
		this->Snow = std::vector<double>(end, 0.0);

#pragma omp parallel for num_threads(this->ParallelNumThread)
		for(signed long nc=0; nc<NSubConfig; nc++)
		{
			for(signed long i=0; i<end; i++)
			{
				double tr=0, ti=0;//real and imag parts of result
				for(size_t j=nc*Np; j<(nc+1)*Np; j++)
				{
					GeometryVector c=pConfig->GetCartesianCoordinates(j);
					tr+=std::cos(this->constraints[i].k.Dot(c));
					ti+=std::sin(this->constraints[i].k.Dot(c));
				}
				this->RhoReal[nc*end+i]=tr;
				this->RhoImag[nc*end+i]=ti;

#pragma omp atomic
				this->Snow[i] += (tr*tr+ti*ti);
			}
		}
		for (signed long i = 0; i < end; i++)
			this->Snow[i] /= (Np*NSubConfig);
		this->RhoValid = true;
	}
	virtual void GetRho(void)
	{
		GetRho_MultiConfig();
	}
	virtual double Energy() //return the total potential energy of the given configuration
	{
		if(this->RhoValid==false)
			this->GetRho_MultiConfig();

		double result=0;
		std::vector<double> threadresult(this->ParallelNumThread, 0.0);

		signed long end=this->constraints.size();

#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long i=0; i<end; i++)
			{
				double s = this->Snow[i];
				double deltaS=(s-this->S0[i]);
				threadresult[numthread]+=this->constraints[i].V*deltaS*deltaS;
			}
		}
		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result+=*iter;

		result+=this->Shift;

		return result;
	}
	virtual void Force(GeometryVector & result, size_t i) //calculate the force for the ith Particle
	{
		if(this->RhoValid==false)
			this->GetRho_MultiConfig();

		result=GeometryVector(this->Dimension);
		std::vector<GeometryVector> threadresult(this->ParallelNumThread, result);

		signed long end=this->constraints.size();

		size_t inc=i/this->Np;//the number of subconfig that contains particle i
#pragma omp parallel num_threads(this->ParallelNumThread) default(shared)
		{
			int numthread=omp_get_thread_num();
#pragma omp for schedule(static)
			for(signed long j=0; j<end; j++)
			{
				GeometryVector k=this->constraints[j].k;
				GeometryVector r=this->pConfig->GetCartesianCoordinates(i);
				double s = this->Snow[j];
				double deltaS=(s-this->S0[j]);
				threadresult[numthread].MinusFrom((4.0/this->Np)*this->constraints[j].V*deltaS*(this->RhoImag[j+inc*end]*std::cos(k.Dot(r))-this->RhoReal[j+inc*end]*std::sin(k.Dot(r)))*k);
			}
		}

		for(auto iter=threadresult.begin(); iter!=threadresult.end(); iter++)
			result.AddFrom(*iter);
		result.MultiplyFrom(1.0/this->NSubConfig);
	}
	virtual void SetConfiguration(const Configuration & Position)
	{
		assert(this->S0.size()==this->constraints.size());
		size_t N=Position.NumParticle();
		if(N%NSubConfig!=0)
		{
			std::cerr<<"Error in ShiftedCCPotential_S2_MultiConfig : the total number of particles is not multiple of the number of configurations!\n";
			std::cerr<<"N="<<N<<", NSubConfig="<<NSubConfig<<'\n';
			assert(false);
		}
		this->Np=N/this->NSubConfig;
		this->ShiftedCCPotential::SetConfiguration(Position);
		this->Volume=this->pConfig->PeriodicVolume();
	}
};


//convert a collective-coordinate potential to a form suitable for Monte Carlo simulations
//require that Energy() only depends on rho (and not depend on r explicitly)
template<typename T>
class ShiftedCCPotential_ForMCMixin : public T, public forMCPotential
{
protected:
	std::vector<double> RhoReal_AfterMove, RhoImag_AfterMove;

public:
	ShiftedCCPotential_ForMCMixin(const T & src) :T(src)
	{
	}

	//parameters : Cartesian coordinates before and after move
	//the potential must have been set to the correct configuration before this function is called
	//return : delta energy
	virtual double TryMove(GeometryVector & prevCart, GeometryVector & afterCart)
	{
		this->RhoReal_AfterMove.resize(this->constraints.size());
		this->RhoImag_AfterMove.resize(this->constraints.size());

		double PrevE = this->Energy();
		signed long end=this->constraints.size();
#pragma omp parallel for num_threads(this->ParallelNumThread)
		for(signed long i=0; i<end; i++)
		{
			this->RhoReal_AfterMove[i]=this->RhoReal[i]-std::cos(this->constraints[i].k.Dot(prevCart))+std::cos(this->constraints[i].k.Dot(afterCart));
			this->RhoImag_AfterMove[i]=this->RhoImag[i]-std::sin(this->constraints[i].k.Dot(prevCart))+std::sin(this->constraints[i].k.Dot(afterCart));
		}
		std::swap(this->RhoImag, this->RhoImag_AfterMove);
		std::swap(this->RhoReal, this->RhoReal_AfterMove);
		double AfterE = this->Energy();
		std::swap(this->RhoImag, this->RhoImag_AfterMove);
		std::swap(this->RhoReal, this->RhoReal_AfterMove);
		return AfterE - PrevE;
	}
	virtual void AcceptMove(void)
	{
		this->RhoImag=this->RhoImag_AfterMove;
		this->RhoReal=this->RhoReal_AfterMove;
	}
};

typedef ShiftedCCPotential_ForMCMixin<ShiftedCCPotential> ShiftedCCPotential_ForMC;

#include <gsl/gsl_linalg.h>
//this class allows box deformation for ShiftedCCPotential or its derivative classes, provides:
//k point list management
//user-specified tildeV(\mathbf k)
//d \mathbf b_i/da_{lm} and d \mathbf k_i/da_{lm}, where b_i is the ith reciprocal basis vector, k_i is the ith k point, 
// and a_{lm} is BasisVector[l].x[m]
//
//does NOT provide:
//management of ShiftedCCPotential::Shift
//implementation of Potential::EnergyDerivativeToBasisVectors 
template<typename Base>
class ShiftedCCPotential_varBoxable : public Base
{
protected:
	double rK;
	GeometryVector CachedBasisVector[::MaxDimension];
	GeometryVector CachedReciprocalBasisVector[::MaxDimension];
	bool SameBasisVector;
	//a k point is made of c[i]*b[i]
	struct kCoeff
	{
		double c[::MaxDimension];
		kCoeff()
		{
			for (int i = 0; i < ::MaxDimension; i++)
				c[i] = 0;
		}
	};
	std::vector<kCoeff> kCoeffs;
	virtual double tildeVk(const GeometryVector & k) = 0;//Eq. (2) in paper-366
	virtual GeometryVector NablatildeVk(const GeometryVector & k) = 0;//gradient of tildeVk
	std::vector<GeometryVector> dBi(DimensionType l, DimensionType m)//d \mathbf b_i / da_{ lm }
	{
		std::vector<GeometryVector> result;
		for (DimensionType i = 0; i < this->Dimension; i++)
		{
			gsl_matrix * a = gsl_matrix_alloc(this->Dimension, this->Dimension);
			for (DimensionType i = 0; i < this->Dimension; i++)
				for (DimensionType j = 0; j < this->Dimension; j++)
					gsl_matrix_set(a, i, j, CachedBasisVector[i].x[j]);
			gsl_vector * x = gsl_vector_alloc(this->Dimension);
			gsl_vector * b = gsl_vector_calloc(this->Dimension);
			gsl_vector_set(b, l, (-1.0)*CachedReciprocalBasisVector[i].x[m]);
			gsl_permutation * p = gsl_permutation_alloc(this->Dimension);
			int s;
			gsl_linalg_LU_decomp(a, p, &s);
			gsl_linalg_LU_solve(a, p, b, x);
			GeometryVector res(this->Dimension);
			for (DimensionType i = 0; i < this->Dimension; i++)
				res.x[i] = gsl_vector_get(x, i);
			result.push_back(res);

			gsl_permutation_free(p);
			gsl_vector_free(b);
			gsl_vector_free(x);
			gsl_matrix_free(a);
		}

		return result;
	}
	std::vector<GeometryVector> dki(DimensionType l, DimensionType m)//d \mathbf k_i / da_{ lm }
	{
		std::vector<GeometryVector> result;
		std::vector<GeometryVector> db = dBi(l, m);
		for (size_t i = 0; i < this->constraints.size(); i++)
		{
			GeometryVector res(this->Dimension);
			for (DimensionType j = 0; j < this->Dimension; j++)
				res.AddFrom(kCoeffs[i].c[j] * db[j]);
			result.push_back(res);
		}
		return result;
	}
public:
	ShiftedCCPotential_varBoxable(DimensionType d, double rK) : Base(d), rK(rK)
	{
		this->SameBasisVector = false;
		this->Base::Shift = 0.0;
		for (DimensionType i = 0; i < MaxDimension; i++)
			CachedBasisVector[i] = GeometryVector(d);
	}
	virtual void SetConfiguration(const Configuration & config)
	{
		//this function does not set Base::Shift
		DimensionType d = this->Potential::Dimension;
		//test if basis vector has never changed
		//if so, we don't need to re-generate the constrained k points
		for (DimensionType i = 0; i < d; i++)
		{
			if (this->CachedBasisVector[i] != config.GetBasisVector(i))
			{
				SameBasisVector = false;
				break;
			}
		}
		if (SameBasisVector)
		{
			this->Base::SetConfiguration(config);
			return;
		}

		GeometryVector b[::MaxDimension];
		for (DimensionType i = 0; i<d; i++)
			b[i] = config.GetReciprocalBasisVector(i);
		PeriodicCellList<Empty> RecipLattice(d, b, 1e100, true);
		RecipLattice.Insert(Empty(), SameDimensionZeroVector(b[0]));
		std::vector<GeometryVector> KPoints;
		std::vector<kCoeff> & temp = this->kCoeffs;
		temp.clear();
		auto IterateFunction = [&KPoints, &d, &temp](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
		{
			bool NonTrivial = false;//is not vector zero and is not symmetrical to another vector
			for (DimensionType i = 0; i<d; i++)
			{
				if (PeriodicShift[i]>0)
				{
					NonTrivial = true;
					break;
				}
				else if (PeriodicShift[i]<0)
				{
					NonTrivial = false;
					break;
				}
			}
			if (NonTrivial)
			{
				KPoints.push_back(shift);
				kCoeff a;
				for (DimensionType i = 0; i < d; i++)
					a.c[i] = PeriodicShift[i];
				temp.push_back(a);
			}
		};
		//Must not refine basis vectors. Doing so will mess up correspondence between constraints and kCoeff
		//RecipLattice.TryRefineBasisVectors();
		RecipLattice.IterateThroughNeighbors(SameDimensionZeroVector(b[0]), rK, IterateFunction);

		this->Base::constraints.clear();
		for (auto iter = KPoints.begin(); iter != KPoints.end(); iter++)
		{
			double Vk = this->tildeVk(*iter);
			this->Base::constraints.push_back(ShiftedCCPotential::KPointConstraint(*iter, Vk));
		}
		this->Base::SetConfiguration(config);

		for (DimensionType i = 0; i < d; i++)
		{
			this->CachedBasisVector[i] = config.GetBasisVector(i);
			this->CachedReciprocalBasisVector[i] = config.GetReciprocalBasisVector(i);
		}
		//carried from old code, however, here one cannot sort them without changing order of kCoeffs
		//auto CompareFunction = [](const GeometryVector & left, const GeometryVector & right) ->bool
		//{
		//	return left.Modulus2()<right.Modulus2();
		//};
		//std::sort(KPoints.begin(), KPoints.end(), CompareFunction);
	}
};

//Potential in Eq. (6) of paper-366, need implementation of:
//protected:
//virtual double tildeVk(const GeometryVector & k);//Eq. (2) in paper-366
//virtual GeometryVector NablatildeVk(const GeometryVector & k);//gradient of tildeVk
class ShiftedCCPotential_varBox : public ShiftedCCPotential_varBoxable<ShiftedCCPotential>
{
protected:
	size_t nbr;
public:
	ShiftedCCPotential_varBox(DimensionType d, double rK = 1.0) : ShiftedCCPotential_varBoxable<ShiftedCCPotential>(d, rK)
	{}
	virtual void SetConfiguration(const Configuration & config)
	{
		ShiftedCCPotential_varBoxable<ShiftedCCPotential>::SetConfiguration(config);
		this->GetRho();

		if (!SameBasisVector || nbr != config.NumParticle())
		{
			//we need to re-calculate this->Shift of basis vector of nbr has changed
			//calculate SumVk
			double SumVk = 0.0;
			for (auto iter = this->constraints.begin(); iter != this->constraints.end(); ++iter)
				SumVk += iter->V;

			nbr = config.NumParticle();
			double v = config.PeriodicVolume();

			this->Shift = nbr*(nbr - 1) / v / 2.0 - nbr / v*SumVk;
		}
	}
	virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	{
		double vF = this->Volume;
		for (DimensionType n = 0; n<Dimension; n++)
			for (DimensionType m = 0; m < Dimension; m++)
			{
				std::vector<GeometryVector> dk = this->dki(n, m);
				double temp = 0.0;
				for (size_t i = 0; i < this->constraints.size(); i++)
					temp += this->NablatildeVk(this->constraints[i].k).Dot(dk[i])*(RhoReal[i] * RhoReal[i] + RhoImag[i] * RhoImag[i] - nbr);
				temp /= vF;

				grad[n*Dimension + m] = temp;
			}
		//get the inverse basis vector matrix
		double phi = this->Energy();
		gsl_matrix * minvbasis = nullptr;
		gsl_matrix * mbasis = gsl_matrix_alloc(Dimension, Dimension);
		minvbasis = gsl_matrix_alloc(Dimension, Dimension);
		for (DimensionType i = 0; i<Dimension; i++)
		{
			GeometryVector temp = this->pConfig->GetBasisVector(i);
			for (DimensionType j = 0; j<Dimension; j++)
				mbasis->data[i*Dimension + j] = temp.x[j];
		}
		gsl_permutation * p = gsl_permutation_calloc(Dimension);
		int signum = 0;
		gsl_linalg_LU_decomp(mbasis, p, &signum);
		gsl_linalg_LU_invert(mbasis, p, minvbasis);
		gsl_permutation_free(p);
		gsl_matrix_free(mbasis);
		for (size_t j = 0; j<Dimension; j++)
		{
			for (size_t k = 0; k<Dimension; k++)
			{
				grad[j*Dimension + k] += (Pressure - phi / vF)*Volume*gsl_matrix_get(minvbasis, k, j);
			}
		}
		gsl_matrix_free(minvbasis);
	}
};

template<typename T>
class K2Template : public T
{
protected:
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1.0 - kl;
		return temp*temp;
	}
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1.0 - kl;
		return ((-2.0)*temp / kl)*k;
	}
public:
	K2Template(DimensionType d) : T(d)
	{}
};

//V(k)=(1/k-1)^2
template<typename T> class KM2Template : public T
{
protected:
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return temp*temp;
	}
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return (-2.0)*temp / kl / kl / kl*k;
	}
public:
	KM2Template(DimensionType d) : T(d)
	{}
};
//V(k)=(1/k-1)^3
template<typename T> class KM3Template : public T
{
protected:
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return temp*temp*temp;
	}
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return (-3.0)*temp*temp / kl / kl / kl*k;
	}
public:
	KM3Template(DimensionType d) : T(d)
	{}
};
//V(k)=(1/k-1)^4
template<typename T> class KM4Template : public T
{
protected:
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return temp*temp*temp*temp;
	}
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return (-4.0)*temp*temp*temp / kl / kl / kl*k;
	}
public:
	KM4Template(DimensionType d) : T(d)
	{}
}; 

template<typename T>
class LemniscateVkTemplate : public T
{
protected:
	double Rc2;
	virtual double tildeVk(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double k2 = k.Modulus2();
		double x2my2 = k.x[0] * k.x[0] - k.x[1] * k.x[1];
		double ratio = k2*k2 / Rc2 / x2my2;
		if (ratio < 1 && x2my2>0)
		{
			double temp = ratio - 1;
			//return temp*temp;
			return temp*temp*temp*temp;
		}
		else
			return 0.0;
	}
	virtual GeometryVector NablatildeVk(const GeometryVector & k)//gradient of tildeVk
	{
		double k2 = k.Modulus2();
		double x2my2 = k.x[0] * k.x[0] - k.x[1] * k.x[1];
		double ratio = k2*k2 / Rc2 / x2my2;
		if (ratio < 1 && x2my2>0)
		{
			double temp = ratio - 1;
			double commonpart = 2 * temp / Rc2 / x2my2 / x2my2;
			//return GeometryVector(4 * k.x[0] * k2*x2my2 - 2 * k.x[0] * k2*k2, 4 * k.x[1] * k2*x2my2 + 2 * k.x[1] * k2*k2)*commonpart;
			return GeometryVector(4 * k.x[0] * k2*x2my2 - 2 * k.x[0] * k2*k2, 4 * k.x[1] * k2*x2my2 + 2 * k.x[1] * k2*k2)*(commonpart*2*temp*temp);
		}
		else
			return GeometryVector(0.0, 0.0);
	}
public:
	LemniscateVkTemplate(double Rc) : T(2, Rc), Rc2(Rc*Rc)
	{}
};



typedef K2Template<ShiftedCCPotential_varBox> K2Potential;
typedef KM3Template<ShiftedCCPotential_varBox> KM3Potential;

//// \phi = \sum_{\mathbf k} tildeV(\mathbf k) ( | \rho(\mathbf k) |^2/N - S0(\mathbf k) )^2, need implementation of:
//protected:
//virtual double tildeVk(const GeometryVector & k);
//virtual GeometryVector NablatildeVk(const GeometryVector & k);//gradient of tildeVk
//virtual double S0k(const GeometryVector & k);
//virtual GeometryVector NablaS0k(const GeometryVector & k);//gradient of S0k
class ShiftedCCPotential_S2_varBox : public ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>
{
protected:
	virtual double S0k(const GeometryVector & k) = 0;
	virtual GeometryVector NablaS0k(const GeometryVector & k) = 0;//gradient of S0k
public:
	ShiftedCCPotential_S2_varBox(DimensionType d, double rK = 1.0) : ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>(d, rK)
	{}
	virtual void SetConfiguration(const Configuration & config)
	{
		ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>::SetConfiguration(config);
		this->GetRho();

		if (!SameBasisVector)
		{
			//we need to re-calculate S0 of basis vector has changed
			this->S0.resize(this->constraints.size());
			for (size_t i = 0; i < this->constraints.size(); i++)
				this->S0[i] = S0k(this->constraints[i].k);
		}
	}
	virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	{
		for (DimensionType n = 0; n < Dimension; n++)
			for (DimensionType m = 0; m < Dimension; m++)
			{
				std::vector<GeometryVector> dk = this->dki(n, m);
				double temp = 0.0;
				for (size_t i = 0; i < this->constraints.size(); i++)
				{
					double ds = ((RhoReal[i] * RhoReal[i] + RhoImag[i] * RhoImag[i]) / N - S0[i]);
					temp += this->NablatildeVk(this->constraints[i].k).Dot(dk[i])*ds*ds;
					temp += this->constraints[i].V*(-2.0)*ds*NablaS0k(this->constraints[i].k).Dot(dk[i]);
				}
				grad[n*Dimension + m] = temp;
			}
		//get the inverse basis vector matrix
		//double phi = this->Energy();
		gsl_matrix * minvbasis = nullptr;
		gsl_matrix * mbasis = gsl_matrix_alloc(Dimension, Dimension);
		minvbasis = gsl_matrix_alloc(Dimension, Dimension);
		for (DimensionType i = 0; i<Dimension; i++)
		{
			GeometryVector temp = this->pConfig->GetBasisVector(i);
			for (DimensionType j = 0; j<Dimension; j++)
				mbasis->data[i*Dimension + j] = temp.x[j];
		}
		gsl_permutation * p = gsl_permutation_calloc(Dimension);
		int signum = 0;
		gsl_linalg_LU_decomp(mbasis, p, &signum);
		gsl_linalg_LU_invert(mbasis, p, minvbasis);
		gsl_permutation_free(p);
		gsl_matrix_free(mbasis);
		for (size_t j = 0; j<Dimension; j++)
		{
			for (size_t k = 0; k<Dimension; k++)
			{
				grad[j*Dimension + k] += (Pressure)*Volume*gsl_matrix_get(minvbasis, k, j);
			}
		}
		gsl_matrix_free(minvbasis);
	}
};

//// \phi = \sum_{\mathbf k} tildeV(\mathbf k) ( | \rho(\mathbf k) |^2/N - S0(\mathbf k) )^2-V0, 
//where V0 is given by Eq. (28) of paper-309.
//need implementation of:
//protected:
//virtual double tildeVk(const GeometryVector & k);
//virtual GeometryVector NablatildeVk(const GeometryVector & k);//gradient of tildeVk
//virtual double S0k(const GeometryVector & k);
//virtual GeometryVector NablaS0k(const GeometryVector & k);//gradient of S0k
class ShiftedCCPotential_S2_varBox_noV0 : public ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>
{
protected:
	virtual double S0k(const GeometryVector & k) = 0;
	virtual GeometryVector NablaS0k(const GeometryVector & k) = 0;//gradient of S0k
public:
	ShiftedCCPotential_S2_varBox_noV0(DimensionType d, double rK = 1.0) : ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>(d, rK)
	{}
	virtual void SetConfiguration(const Configuration & config)
	{
		ShiftedCCPotential_varBoxable<ShiftedCCPotential_S2>::SetConfiguration(config);
		this->GetRho();

		if (!SameBasisVector)
		{
			//we need to re-calculate S0 of basis vector has changed
			this->S0.resize(this->constraints.size());
			this->Shift = 0.0;
			for (size_t i = 0; i < this->constraints.size(); i++)
			{
				double s0 = S0k(this->constraints[i].k);
				this->S0[i] = s0;
				this->Shift -= (s0 - 1)*(s0 - 1)*this->constraints[i].V;
			}
		}
	}
	virtual void EnergyDerivativeToBasisVectors(double * grad, double Pressure)
	{
		for (DimensionType n = 0; n < Dimension; n++)
			for (DimensionType m = 0; m < Dimension; m++)
			{
				std::vector<GeometryVector> dk = this->dki(n, m);
				double temp = 0.0;
				for (size_t i = 0; i < this->constraints.size(); i++)
				{
					double ds = ((RhoReal[i] * RhoReal[i] + RhoImag[i] * RhoImag[i]) / N - S0[i]);
					temp += this->NablatildeVk(this->constraints[i].k).Dot(dk[i])*(ds*ds-(S0[i]-1.0)*(S0[i]-1.0)) + this->constraints[i].V*(1.0-S0[i] - ds)*2.0*NablaS0k(this->constraints[i].k).Dot(dk[i]);
				}
				grad[n*Dimension + m] = temp;
			}
		//get the inverse basis vector matrix
		//double phi = this->Energy();
		gsl_matrix * minvbasis = nullptr;
		gsl_matrix * mbasis = gsl_matrix_alloc(Dimension, Dimension);
		minvbasis = gsl_matrix_alloc(Dimension, Dimension);
		for (DimensionType i = 0; i<Dimension; i++)
		{
			GeometryVector temp = this->pConfig->GetBasisVector(i);
			for (DimensionType j = 0; j<Dimension; j++)
				mbasis->data[i*Dimension + j] = temp.x[j];
		}
		gsl_permutation * p = gsl_permutation_calloc(Dimension);
		int signum = 0;
		gsl_linalg_LU_decomp(mbasis, p, &signum);
		gsl_linalg_LU_invert(mbasis, p, minvbasis);
		gsl_permutation_free(p);
		gsl_matrix_free(mbasis);
		for (size_t j = 0; j<Dimension; j++)
		{
			for (size_t k = 0; k<Dimension; k++)
			{
				grad[j*Dimension + k] += (Pressure)*Volume*gsl_matrix_get(minvbasis, k, j);
			}
		}
		gsl_matrix_free(minvbasis);
	}
};


//S0(k)=k^2
template<typename T> class S0K2Template : public T
{
protected:
	virtual double S0k(const GeometryVector & k)
	{
		return k.Modulus2();
	}
	virtual GeometryVector NablaS0k(const GeometryVector & k)
	{
		return 2.0*k;
	}
public:
	S0K2Template(DimensionType d) : T(d)
	{}
};

//S0(k)=k^n
template<typename T> class S0KnTemplate : public T
{
protected:
	double n;
	virtual double S0k(const GeometryVector & k)
	{
		return std::pow(k.Modulus2(), 0.5*n);
	}
	virtual GeometryVector NablaS0k(const GeometryVector & k)
	{
		return n*std::pow(k.Modulus2(), 0.5*n-1)*k;
	}
public:
	S0KnTemplate(DimensionType d, double n) : T(d), n(n)
	{}
};

//S0(k)=(k^n+1)/2
template<typename T> class S0KlcTemplate : public T
{
protected:
	double n;
	virtual double S0k(const GeometryVector & k)
	{
		return 0.5*(std::pow(k.Modulus2(), 0.5*n)+1);
	}
	virtual GeometryVector NablaS0k(const GeometryVector & k)
	{
		return 0.5*n*std::pow(k.Modulus2(), 0.5*n - 1)*k;
	}
public:
	S0KlcTemplate(DimensionType d, double n) : T(d), n(n)
	{}
};

template<typename T> class S0KM2Template : public T
{
protected:
	virtual double S0k(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return temp*temp;
	}
	virtual GeometryVector NablaS0k(const GeometryVector & k)//gradient of tildeVk
	{
		double kl = std::sqrt(k.Modulus2());
		double temp = 1 / kl - 1;
		return (-2.0)*temp / kl / kl / kl*k;
	}
public:
	S0KM2Template(DimensionType d) : T(d)
	{}
};


template<typename T> class S0KMLogTemplate : public T
{
protected:
	virtual double S0k(const GeometryVector & k)//Eq. (2) in paper-366
	{
		double kl2 = k.Modulus2();
		return -0.5*std::log(kl2);
	}
	virtual GeometryVector NablaS0k(const GeometryVector & k)//gradient of tildeVk
	{
		return k*(-1.0/k.Modulus2());
	}
public:
	S0KMLogTemplate(DimensionType d) : T(d)
	{}
};



typedef S0K2Template<K2Template<ShiftedCCPotential_S2_varBox> > K2K2Potential;
typedef S0K2Template<K2Template<ShiftedCCPotential_S2_varBox_noV0> > K2K2Potential_noV0;
typedef S0K2Template<KM3Template<ShiftedCCPotential_S2_varBox_noV0> > KM3K2Potential_noV0;


typedef S0K2Template<KM2Template<ShiftedCCPotential_S2_varBox> > KM2K2Potential;
typedef S0K2Template<KM3Template<ShiftedCCPotential_S2_varBox> > KM3K2Potential;
typedef S0K2Template<KM4Template<ShiftedCCPotential_S2_varBox> > KM4K2Potential;

typedef S0KM2Template<KM3Template<ShiftedCCPotential_S2_varBox> > KM3KM2Potential;

typedef S0KMLogTemplate<KM3Template<ShiftedCCPotential_S2_varBox> > KM3KMLogPotential;

typedef S0KnTemplate<KM2Template<ShiftedCCPotential_S2_varBox> > KM2KnPotential;
typedef S0KnTemplate<KM3Template<ShiftedCCPotential_S2_varBox> > KM3KnPotential;
typedef S0KnTemplate<KM4Template<ShiftedCCPotential_S2_varBox> > KM4KnPotential;

typedef S0KlcTemplate<KM2Template<ShiftedCCPotential_S2_varBox> > KM2KlcPotential;
typedef S0KlcTemplate<KM3Template<ShiftedCCPotential_S2_varBox> > KM3KlcPotential;
typedef S0KlcTemplate<KM4Template<ShiftedCCPotential_S2_varBox> > KM4KlcPotential;

#endif
