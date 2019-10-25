#include "Potential.h"
#include <cmath>
#include <cassert>

const double Potential_Check_Delta = 1e-7;
void Potential::Check(Configuration c, size_t Nbr)
{
	this->SetConfiguration(c);
	GeometryVector temp;
	GeometryVector force, force2, force3;
	this->Force(force, Nbr);

	std::vector<GeometryVector> fs;
	this->AllForce(fs);
	force2 = fs[Nbr];
	double E0Other=this->AllForceAndEnergy(fs);
	force3 = fs[Nbr];

	double E0 = this->Energy();
	std::cout << "Energy() vs. AllForceAndEnergy() check : " << E0 / E0Other << std::endl;
	temp = c.GetCartesianCoordinates(Nbr);
	temp.x[0] += Potential_Check_Delta;
	c.MoveParticle(Nbr, c.CartesianCoord2RelativeCoord(temp));
	this->SetConfiguration(c);
	double E1 = this->Energy();
	std::cout << "force check : " << (E0 - E1) / Potential_Check_Delta / force.x[0] << '\n';
	std::cout << "AllForce check : " << (E0 - E1) / Potential_Check_Delta / force2.x[0] << '\n';
	std::cout << "AllForceAndEnergy check : " << (E0 - E1) / Potential_Check_Delta / force3.x[0] << '\n';
}
//this function will print a number very close to 1 when EnergyDerivativeToBasisVectors() is consistent with Energy(), and c is not a local minimum of this potential
void Potential::Check2(Configuration c, size_t Nbr, double Pressure)
{
	size_t d2 = Dimension*Dimension;
	Nbr = Nbr%d2;
	double * der = new double[d2];
	this->SetConfiguration(c);
	this->EnergyDerivativeToBasisVectors(der, Pressure);
	double E0 = this->Energy() + Pressure*c.PeriodicVolume();

	GeometryVector bas[::MaxDimension];
	for (DimensionType i = 0; i < Dimension; i++)
		bas[i] = c.GetBasisVector(i);
	bas[Nbr / Dimension].x[Nbr%Dimension] += Potential_Check_Delta;
	c.ChangeBasisVector(bas);
	this->SetConfiguration(c);
	double E1 = this->Energy() + Pressure*c.PeriodicVolume();
	std::cout << (E1 - E0) / Potential_Check_Delta / der[Nbr] << '\n';

	delete[] der;
}
Potential::~Potential()
{
}
PairPotential::PairPotential(DimensionType dimension, double Rcut):Potential(dimension)
{
	if(dimension<=0) throw "error in PairPotential::PairPotential : unsupported dimension";
	if (Rcut < 0)
	{
		std::cerr<< "error in PairPotential::PairPotential : unsupported Rcut="<<Rcut<<'\n';
	}

	this->Dimension=dimension;
	this->Rcut2=Rcut*Rcut;
	this->Rcut=Rcut;
	this->pConfig=nullptr;
}
PairPotential::~PairPotential()
{
}
void PairPotential::SetConfiguration(const Configuration & Position)
{
	this->pConfig = & Position;
}
double PairPotential::Energy()
{
	size_t Np = this->pConfig->NumParticle();
	//pre-allocation of internal neighbor-cell cache in Configuration
	this->pConfig->PrepareIterateThroughNeighbors(this->Rcut);

	assert(this->pConfig != nullptr);
	const Configuration & conf = *this->pConfig;
	double Energy = 0;
#pragma omp parallel num_threads(this->ParallelNumThread)
	{
		KahanAccumulation threadSum;
#pragma omp for 
		for (long i = 0; i < conf.NumParticle(); i++)
		{
			PairPotential & backup = (*this);
			AtomInfo a1 = pConfig->GetCharacteristics(i);
			const Configuration * pConfig2 = this->pConfig;
			//Configuration::particle * trial = this->pConfig->GetParticle(i);
			conf.IterateThroughNeighbors(i, this->Rcut, [&threadSum, &backup, &i, &a1, &pConfig2](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void
			{
				if (SourceAtom != i || LatticeShift.Modulus2() != 0)
				{
					AtomInfo a2 = pConfig2->GetCharacteristics(SourceAtom);
					double temp = backup.EnergyFormula(shift, a1, a2);
					//#pragma omp atomic
					//					Energy += temp;
					threadSum = KahanSum(threadSum, temp);
				}
			});
		}
#pragma omp atomic
		Energy += threadSum.sum;
	}
	Energy /= 2;
	return Energy;
}

double PairPotential::Energy(const Configuration & conf)
{
	double Energy=0;
	for(Configuration::Index i=0; i<conf.NumParticle(); i++)
	{
		AtomInfo a1=conf.GetCharacteristics(i);
		PairPotential & backup = ( *this);
		const Configuration * pConfig2 = &conf;
		conf.IterateThroughNeighbors(i, this->Rcut, [&Energy, &backup, &i, &a1, &pConfig2](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void 
		{
			if(SourceAtom!=i || LatticeShift.Modulus2()!=0) 
			{
				AtomInfo a2=pConfig2->GetCharacteristics(SourceAtom);
				double temp = backup.EnergyFormula(shift, a1, a2);
				Energy+=temp;
			}
		});
	}
	Energy/=2;
	return Energy;
}
void PairPotential::Force(GeometryVector & result, size_t i)
{
	assert(this->pConfig!=nullptr);
	result=GeometryVector(this->Dimension);

	PairPotential & backup = ( *this);
	//Configuration::particle * trial = this->pConfig->GetParticle(i);
		AtomInfo a1=pConfig->GetCharacteristics(i);
		const Configuration * pConfig2=this->pConfig;
	this->pConfig->IterateThroughNeighbors(i, this->Rcut, [&result, &backup, &i, &a1, &pConfig2](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void 
	{
		if(SourceAtom==i)
			return;
		GeometryVector temp=SameDimensionZeroVector(shift);
				AtomInfo a2=pConfig2->GetCharacteristics(SourceAtom);

		//if(shift.Modulus2()> ::LengthPrecision* ::LengthPrecision) 
			backup.ForceFormula((-1.0)*shift, temp, a1, a2);
		result.AddFrom(temp);
	}
	);
	return;
}
void PairPotential::Force(const Configuration & conf, GeometryVector & result, size_t i)
{
	result = GeometryVector(this->Dimension);

	PairPotential & backup = (*this);
	AtomInfo a1 = conf.GetCharacteristics(i);
	const Configuration * pConfig2 = &conf;
	conf.IterateThroughNeighbors(i, this->Rcut, [&result, &backup, &i, &a1, &pConfig2](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void
	{
		if (SourceAtom == i)
			return;
		GeometryVector temp = SameDimensionZeroVector(shift);
		AtomInfo a2 = pConfig2->GetCharacteristics(SourceAtom);

		//if(shift.Modulus2()> ::LengthPrecision* ::LengthPrecision) 
		backup.ForceFormula((-1.0)*shift, temp, a1, a2);
		result.AddFrom(temp);
	}
	);
	return;
}
void PairPotential::AllForce(std::vector<GeometryVector> & results)
{
	size_t Np=this->pConfig->NumParticle();
	results.resize(Np);
	//pre-allocation of internal neighbor-cell cache in Configuration
	this->pConfig->PrepareIterateThroughNeighbors(this->Rcut);
#pragma omp parallel for num_threads(this->ParallelNumThread)
	for(long i=0; i<Np; i++)
		this->Force(results[i], i);
}

double PairPotential::ForceAndEnergyFormula(const GeometryVector & offset, GeometryVector & result, const AtomInfo & a1, const AtomInfo & a2)
{
	//default implementation
	this->ForceFormula(offset, result, a1, a2);
	return this->EnergyFormula(offset, a1, a2);
}

double PairPotential::AllForceAndEnergy(std::vector<GeometryVector> & results)
{
	size_t Np = this->pConfig->NumParticle();
	results = std::vector<GeometryVector>(Np, GeometryVector(this->Dimension));
	//pre-allocation of internal neighbor-cell cache in Configuration
	this->pConfig->PrepareIterateThroughNeighbors(this->Rcut);

	assert(this->pConfig != nullptr);
	const Configuration & conf = *this->pConfig;
	double Energy = 0;
#pragma omp parallel num_threads(this->ParallelNumThread)
	{
		KahanAccumulation threadSum;
#pragma omp for schedule(guided)
		for (long i = 0; i < conf.NumParticle(); i++)
		{
			PairPotential & backup = (*this);
			AtomInfo a1 = pConfig->GetCharacteristics(i);
			const Configuration * pConfig2 = this->pConfig;
			//Configuration::particle * trial = this->pConfig->GetParticle(i);
			conf.IterateThroughNeighbors(i, this->Rcut, [&threadSum, &backup, &i, &a1, &pConfig2, &results](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void
			{
				if (SourceAtom != i || LatticeShift.Modulus2() != 0)
				{
					AtomInfo a2 = pConfig2->GetCharacteristics(SourceAtom);
					GeometryVector temp2 = SameDimensionZeroVector(shift);
					double temp = backup.ForceAndEnergyFormula((-1.0)*shift, temp2, a1, a2);
					threadSum = KahanSum(threadSum, temp);
					results[i].AddFrom(temp2);
				}
			});
		}
#pragma omp atomic
		Energy += threadSum.sum;
	}
	Energy /= 2;
	return Energy;
}

double PairPotential::SecondDerivative(const GeometryVector & offset, DimensionType i, DimensionType j, const AtomInfo & a1, const AtomInfo & a2)//(d^2V)/(dxi*dxj), recommended to rewrite it
{
	GeometryVector off(offset);
	GeometryVector f1(this->Dimension);
	this->ForceFormula(off, f1, a1, a2);

	off.x[i]+= ::LengthPrecision;
	GeometryVector f2(this->Dimension);
	this->ForceFormula(off, f2, a1, a2);

	return (f1.x[j]-f2.x[j])/ ::LengthPrecision;
}


double IsotropicPairPotential::EnergyFormula(const GeometryVector & offset, const AtomInfo & a1, const AtomInfo & a2)
{
	return this->IsotropicPairEnergy(std::sqrt(offset.Modulus2()), a1, a2);
}

void IsotropicPairPotential::ForceFormula(const GeometryVector & offset, GeometryVector & result, const AtomInfo & a1, const AtomInfo & a2)
{
	double dist=std::sqrt(offset.Modulus2());
	double coeff=this->IsotropicForceFormula(dist, a1, a2)/dist;
	result=coeff*offset;
}

IsotropicPairPotential::IsotropicPairPotential(DimensionType dimension, double Rcut):PairPotential(dimension, Rcut)
{
}
IsotropicPairPotential::~IsotropicPairPotential()
{
}

const double deltaR=1e-10;
double IsotropicPairPotential::IsotropicForceFormula(double distance, const AtomInfo & a1, const AtomInfo & a2)
{
	double e1=this->IsotropicEnergyFormula(distance-deltaR, a1, a2);
	double e2=this->IsotropicEnergyFormula(distance+deltaR, a1, a2);
	return (e1-e2)/(2*deltaR);
}

double IsotropicPairPotential::IsotropicForceAndEnergyFormula(double distance, const AtomInfo & a1, const AtomInfo & a2, double & force)
{
	force = this->IsotropicForceFormula(distance, a1, a2);
	return this->IsotropicEnergyFormula(distance, a1, a2);
}
double IsotropicPairPotential::ForceAndEnergyFormula(const GeometryVector & offset, GeometryVector & result, const AtomInfo & a1, const AtomInfo & a2)
{
	double dist = std::sqrt(offset.Modulus2());
	double force, energy;
	energy = this->IsotropicForceAndEnergyFormula(dist, a1, a2, force);
	double coeff = force / dist;
	result = coeff*offset;
	return energy;
}


double IsotropicPairPotential::SecondDerivative(const GeometryVector & offset, DimensionType i, DimensionType j, const AtomInfo & a1, const AtomInfo & a2)//(d^2V)/(dxi*dxj)
{
	double r2 = offset.Modulus2();
	double r = std::sqrt(r2);
	double coeff2 = offset.x[i] * offset.x[j] / r2;
	double coeff1 = coeff2 / r;
	if (i == j)
		coeff1 -= (1 / r);

	return coeff1*this->IsotropicForceFormula(r, a1, a2) + coeff2*this->IsotropicSecondDerivative(r, a1, a2);
}
double IsotropicPairPotential::IsotropicSecondDerivative(double distance, const AtomInfo & a1, const AtomInfo & a2)//-(dF)/(dx), recommended to override
{
	double e1 = this->IsotropicForceFormula(distance - deltaR, a1, a2);
	double e2 = this->IsotropicForceFormula(distance + deltaR, a1, a2);
	return (e1 - e2) / (2 * deltaR);
}


//Energy derivative respect to basis vectors
//grad[n*dim+m] is the energy derivative of mth (coordinate) component of nth basis vector
void PairPotential::EnergyDerivativeToBasisVectors(double * grad, double Pressure)
{
	assert(this->pConfig!=nullptr);

	DimensionType dim=this->pConfig->GetDimension();
	for(size_t i=0; i<dim*dim; i++)
		grad[i]=0.0;
	double Volume=this->pConfig->PeriodicVolume();
	//pre-allocation of internal neighbor-cell cache in Configuration
	this->pConfig->PrepareIterateThroughNeighbors(this->Rcut);
#pragma omp parallel for num_threads(this->ParallelNumThread)
	for(long i=0; i<this->pConfig->NumParticle(); i++)
	{
		//Configuration::particle * pA=this->pConfig->GetParticle(i);
		const Configuration * pTemp=this->pConfig;
		AtomInfo a1=pConfig->GetCharacteristics(i);
		//derivatives of cell basis vectors
		this->pConfig->IterateThroughNeighbors(i, this->Rcut, 
			[&] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t srcAtom) ->void
		{
			AtomInfo a2=pTemp->GetCharacteristics(srcAtom);
			if(shift.Modulus2()< LengthPrecision*LengthPrecision)
				return;
			GeometryVector RelativeShift = pTemp->GetRelativeCoordinates(srcAtom)-pTemp->GetRelativeCoordinates(i);
			for(DimensionType i=0; i<dim; i++)
				RelativeShift.x[i]+=PeriodicShift[i];
			GeometryVector pairForce(dim);
			this->ForceFormula(shift, pairForce, a1, a2);
			for(DimensionType m=0; m<dim; m++)
				for(DimensionType n=0; n<dim; n++)
				{
					double temp=(-0.5)*(pairForce.x[m])*RelativeShift.x[n];
#pragma omp atomic
					grad[n*dim+m] += temp;
				}
		}
		);
	}
	for(size_t j=0; j<dim; j++)
	{
		for(size_t k=0; k<dim; k++)
		{
			grad[j*dim+k]+=Pressure*Volume*this->pConfig->GetReciprocalBasisVector(j).x[k]/2/pi;
		}
	}
}
void PairPotential::EnergyDerivativeToBasisVectors(const Configuration & conf, double * grad, double Pressure)
{
	this->SetConfiguration(conf);
	return this->EnergyDerivativeToBasisVectors(grad, Pressure);
}

