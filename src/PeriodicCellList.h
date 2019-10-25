#ifndef PeriodicCellList_INCLUDED
#define PeriodicCellList_INCLUDED

#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <list>
#include <stdexcept>

#include <omp.h>

#include "GeometryVector.h"
#include "RandomGenerator.h"
#include "etc.h"

const size_t ParticlePerCell = 1;


struct NeighborListElementType
{
	size_t SourceParticle;
	signed long PeriodicShift[::MaxDimension];
	GeometryVector LatticeShift;
};
template<typename OtherCharacteristics>
class PeriodicCellList
{
	template<typename a> friend class PeriodicCellList;
	//some concepts:
	//there are 3 types of coordinates:
	//1. Cartesian Coordinates
	//2. Relative Coordinates, this is the coordinates relative to the periodic box base vectors.
	//		The periodic box base vectors are PeriodicCellList::BasisVector
	//3. Lattice Coordinates, this is the coordinates relative to the cell lattice(i.e. the base vectors of a cell list cell).
public:
	typedef unsigned long Index;

protected:
	DimensionType Dimension;
	GeometryVector BasisVector[::MaxDimension];//the base vector of a periodic cell
	mutable signed long CellRank[::MaxDimension];//number of cells along each basis vector
	mutable std::vector<std::vector<size_t > > Cells;//content Particles in the cells
	std::vector<GeometryVector> ParticleRelatives;
	std::vector<GeometryVector> ParticleCartesians;
	std::vector<OtherCharacteristics> ParticleCharacteristics;
	double ApproximateCellSize;// the best choice of cell size obtained in constructor when user input cell size and cell rank. We will adjust CellRank so that their size is roughly equal to this quantity.

	//reciprocal basis vector, calculated when needed
	mutable bool ReciprocalBasisVectorValid[::MaxDimension];
	mutable GeometryVector ReciprocalBasisVector[::MaxDimension];

	//member data for vicinity lattice
	struct VicinityLatticesMember
	{
		GeometryVector CartesianCoordinates;
		signed long LatticeShift[MaxDimension];
	};
	mutable std::vector<VicinityLatticesMember> VicinityLattices;
	mutable std::map<double, size_t> VicinityLatticePartitions;
	mutable double VicinityLatticesRc2;
	bool SortedVicinityList;//use sorted vicinity lattice list, which takes some time to sort, but makes iteration more efficient

	void UpdateCartesianCoordinates(size_t a)
	{
		if (!(this->ParticleRelatives[a].x[0] < 1 + ::LengthPrecision) || !(this->ParticleRelatives[a].x[0] >= 0 - ::LengthPrecision))
		{
			std::stringstream ss;
			ss << "Error in PeriodicCellList::UpdateCartesianCoordinates : invalid ParticleRelatives. ParticleRelatives[" << a << "]=" << this->ParticleRelatives[a];
			throw std::invalid_argument(ss.str());
		}
		this->ParticleCartesians[a] = GeometryVector(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->ParticleCartesians[a].AddFrom((this->BasisVector[i])*(this->ParticleRelatives[a].x[i]));
	}
	void PutToPeriodicBox(size_t a)
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->ParticleRelatives[a].x[i] -= std::floor(this->ParticleRelatives[a].x[i]);
	}
	Index ConvertIndex(const signed long * index) const//convert dimensional index to the 1-dimension cell index
	{
		Index result = index[0];
		for (DimensionType i = 1; i<this->Dimension; i++)
		{
			assert(index[i] >= 0);
			assert(index[i]<this->CellRank[i]);
			result *= this->CellRank[i];
			result += index[i];
		}
		return result;
	}
	void GetFullIndex(signed long * index, const GeometryVector & RelativeCoordinates) const
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			index[i] = (signed long)(std::floor(RelativeCoordinates.x[i] * this->CellRank[i]));
			if (index[i]<0 || index[i] >= this->CellRank[i])
			{
				std::stringstream ss;
				ss << "Error in Configuration::GetFullIndex : index invalid. index=" << index[i] << ", Relative Coordinate=" << RelativeCoordinates.x[i] << ", Cell Rank=" << this->CellRank[i] << '\n';
				throw std::invalid_argument(ss.str());
			}
		}
	}
	Index GetIndex(const size_t a) const //get the index of the cell that a should be placed in, the LatticeCoordinates of particle must be in box
	{
		signed long index[::MaxDimension];
		this->GetFullIndex(index, this->ParticleRelatives[a]);
		return this->ConvertIndex(index);
	}
	bool UpdateCellList(std::vector<size_t>::iterator a, std::vector<size_t> & from)//move a to another cell if necessary, from is the cell that a comes from, return true if moved
	{
		//calculate new index
		Index newindex = this->GetIndex((*a));
		if (&(from) == &(this->Cells[newindex]))//don't need to move
			return false;

		//add to new cell
		this->Cells[newindex].push_back(*a);

		//remove from original cell
		*a = from.back();
		from.pop_back();
		return true;
	}
	bool _VicinityLatticeIteration(GeometryVector now, DimensionType DimIDo, GeometryVector nowRelative, GeometryVector * ortho, double Rmax2, double Rmin2) const
	{
		const double VicinityLatticeRmin_Tolerance = 1e-10;
		if (DimIDo == 0)
		{
			double temp = now.Modulus2();
			if (temp<Rmax2)
			{
				bool ToAdd;
				if (temp >(1.0 + VicinityLatticeRmin_Tolerance)*Rmin2)
					ToAdd = true;
				else if (temp < (1.0 - VicinityLatticeRmin_Tolerance)*Rmin2)
					ToAdd = false;
				else
				{
					//cannot tell if the lattice point was already included just from the distance
					//instead, go through the entire VicinityLattices vector to see if it is already added
					VicinityLatticesMember temp;
					for (DimensionType i = 0; i<this->Dimension; i++)
					{
						double & x = nowRelative.x[i];
						temp.LatticeShift[i] = x >= 0 ? static_cast<signed long>(x + 0.5) : static_cast<signed long>(x - 0.5);
					}
					ToAdd = true;
					for (auto iter = this->VicinityLattices.begin(); iter != this->VicinityLattices.end(); iter++)
					{
						bool LatticeShiftIdentical = true;
						for (DimensionType d = 0; d<this->Dimension; d++)
							if (iter->LatticeShift[d] != temp.LatticeShift[d])
							{
								LatticeShiftIdentical = false;
								break;
							}
						if (LatticeShiftIdentical)
						{
							ToAdd = false;
							break;
						}
					}
				}

				if (ToAdd)
				{
					VicinityLatticesMember temp;
					temp.CartesianCoordinates = now;
					for (DimensionType i = 0; i<this->Dimension; i++)
					{
						double & x = nowRelative.x[i];
						temp.LatticeShift[i] = x >= 0 ? static_cast<signed long>(x + 0.5) : static_cast<signed long>(x - 0.5);
					}

					this->VicinityLattices.push_back(temp);
				}
				return true;
			}
			else
				return false;
		}


		GeometryVector k = now;
		for (DimensionType i = 0; i<DimIDo; i++)
		{
			GeometryVector & e = ortho[i];
			k.MinusFrom(e*(k.Dot(e) / e.Dot(e)));
		}
		if (k.Dot(k)>Rmax2)
			return false;
		//right now k is the normal vector of plane (e1, e2,...) in subspace (e1, e2, ... , now)
		GeometryVector a = now - k;//the part of now in that plane
		GeometryVector ethis = this->BasisVector[DimIDo - 1] * (1.0 / this->CellRank[DimIDo - 1]);//base vector that this function will loop through
		for (DimensionType i = 0; i<DimIDo - 1; i++)
		{
			GeometryVector & e = ortho[i];
			a.MinusFrom(e*(a.Dot(e) / e.Dot(e)));
			ethis.MinusFrom(e*(ethis.Dot(e) / e.Dot(e)));
		}
		double coeff = 0;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (std::abs(ethis.x[i])> ::LengthPrecision)
			{
				coeff = a.x[i] / ethis.x[i];
				break;
			}
		}//a equals to ethis*coeff
		assert((a - coeff*ethis).Modulus2()< ::LengthPrecision);

		double CoeffMaxInPlane = std::sqrt((Rmax2 - k.Modulus2()) / ethis.Modulus2());//sqrt(R^2-k^2)/|ethis|
		signed long imin = (signed long)(std::ceil((-1)*CoeffMaxInPlane - coeff));
		signed long imax = (signed long)(std::floor(CoeffMaxInPlane - coeff));

		//int nearest=std::floor(coeff);
		GeometryVector Unit(this->Dimension);
		Unit.x[DimIDo - 1] = 1;
		for (signed long i = imin; i <= imax; i++)
		{
			this->_VicinityLatticeIteration(now + i*this->BasisVector[DimIDo - 1] * (1.0 / this->CellRank[DimIDo - 1]), DimIDo - 1, nowRelative + i*Unit, ortho, Rmax2, Rmin2);
		}
		return true;
	}

	//members related to neighbor list
	mutable bool NeighborListMode;//true: neighbor list mode, false: cell list mode
	mutable std::vector< std::vector<NeighborListElementType> > NeighborList;
	mutable std::vector<GeometryVector> NeighborListStartingCartesians;
	mutable double NeighborListMaxMovement2SinceBuilt;
	mutable double NeighborListRange;

	mutable std::vector<char> NeighborListEnabledForThisParticle;

	//members related to neighbor list with cell change
	mutable GeometryVector NeighborListStartingReciprocalBasis[::MaxDimension];
	mutable double NeighborListCellDeformationCoeff;
	mutable bool NeighborListLatticeShiftValid;

	double getLambdaMaxOfDeformation()
	{
		//find an upper bound for the absolute value of the eigenvalues of B * B_o^{-1}
		//where B is the basis vector, and B_o is NeighborListStartingReciprocalBasis
		//any upper bound gives a safe estimate of NeighborListCellDeformationCoeff
		//we here implement a crute estimate:
		//let BB_o^{-1}=I+delta
		//the largest eigenvalue is bounded by the largest absolute row sum or column sum
		//this gives an upper bound for the eigenvalues of delta
		//add 1, and we have an upper bound for the eigenvalues of BB_o^{-1}
		double max = 0.0;
		for (DimensionType i = 0; i < this->Dimension; i++)
		{
			double sumAbs = 0.0;
			for (DimensionType j = 0; j < this->Dimension; j++)
			{
				double temp = this->BasisVector[i].Dot(this->NeighborListStartingReciprocalBasis[j])/2/pi;
				if (i == j)
					temp -= 1;
				sumAbs += std::abs(temp);
			}
			max = std::max(max, sumAbs);
		}
		return max + 1.0;
	}

	void IterateThroughNeighbors_cellList(const GeometryVector & RelativeCoordinate, double Rc, const std::function<void(const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)> callback, bool * pAbort = nullptr) const
	{
		Rc = std::abs(Rc);
		double Rc2 = Rc*Rc;
		GeometryVector StartRelativeCoordinate = RelativeCoordinate;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			StartRelativeCoordinate.x[i] -= std::floor(StartRelativeCoordinate.x[i]);
			StartRelativeCoordinate.x[i] -= std::floor(StartRelativeCoordinate.x[i]);
		}
		//find the true Rc, which consideres the distance from a particle to a lattice point
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
			else
				maxdist.MinusFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
		}
		double newRc = Rc + std::sqrt(maxdist.Modulus2());
		double newRc2 = newRc*newRc;

		auto IterationEnd = this->VicinityLattices.end();
		//Get new VicinityLattice if necessary
		if (this->SortedVicinityList)
		{
			if ((newRc2)>this->VicinityLatticesRc2)
				this->GetVicinityLattice(newRc*1.1);
			IterationEnd = this->VicinityLattices.end();
		}
		else
		{
			auto ub = this->VicinityLatticePartitions.lower_bound(newRc2);
			auto lb = ub;
			if (this->VicinityLatticePartitions.size()>1)
				lb--;
			if ((newRc2)>this->VicinityLatticesRc2)
			{
				this->GetVicinityLattice(newRc);
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, this->VicinityLattices.size()));
				IterationEnd = this->VicinityLattices.end();
			}
			else if (ub->first == newRc2)
			{
				IterationEnd = ub->second + this->VicinityLattices.begin();
			}
			else if (ub->second - lb->second>10)
			{
				auto titer = std::partition(this->VicinityLattices.begin() + lb->second, this->VicinityLattices.begin() + ub->second,
					[&newRc2](const VicinityLatticesMember & a) ->bool
				{
					return a.CartesianCoordinates.Modulus2()<newRc2;
				});
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, titer - this->VicinityLattices.begin()));
				IterationEnd = titer;
			}
			else
				IterationEnd = ub->second + this->VicinityLattices.begin();
		}

		//Get cartesian coordinate of start point
		GeometryVector StartCartesianCoordinate(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			StartCartesianCoordinate.AddFrom(this->BasisVector[i] * StartRelativeCoordinate.x[i]);

		//Get the (dimensional) index of the cell that the start point is in
		signed long index[::MaxDimension];
		this->GetFullIndex(index, StartRelativeCoordinate);
		//iterate through cells
		for (auto iter = this->VicinityLattices.begin(); iter != IterationEnd; iter++)
		{
			//the vicinity lattice list can be larger than what's needed
			if (this->SortedVicinityList)
			{
				if (iter->CartesianCoordinates.Modulus2()>newRc2)
					break;
			}
			else
			{
				if (iter->CartesianCoordinates.Modulus2()>newRc2)
					continue;
			}


			signed long newindex[::MaxDimension]; //the dimensional index for the cell
			signed long newshift[::MaxDimension]; //the difference of periodic cell of new particle ant current particle, in terms of relative coordinate
			for (DimensionType i = 0; i<this->Dimension; i++)
			{
				signed long n = index[i] + iter->LatticeShift[i];

				if (n >= 0)
				{
					newindex[i] = (n) % this->CellRank[i];
					newshift[i] = (n) / this->CellRank[i];
				}
				else
				{
					newshift[i] = (-1) - (((-1) - n) / this->CellRank[i]);
					newindex[i] = this->CellRank[i] - 1 - (((-1) - n) % this->CellRank[i]);
				}
			}

			const std::vector<size_t > & thiscell = this->Cells[this->ConvertIndex(newindex)];

			if (thiscell.size() != 0)
			{
				GeometryVector LatticeShift(this->Dimension);
				for (int i = 0; i<this->Dimension; i++)
					LatticeShift.AddFrom(static_cast<double>(newshift[i])*this->BasisVector[i]);
				//iterate through Particles
				for (auto iter2 = thiscell.begin(); iter2 != thiscell.end(); iter2++)
				{
					if (pAbort != nullptr)
						if (*pAbort == true)
							return;

					GeometryVector shift = this->ParticleCartesians[*iter2] - StartCartesianCoordinate;

					shift.AddFrom(LatticeShift);

					if (shift.Modulus2()<Rc2)
						callback(shift, LatticeShift, newshift, *iter2);
				}
			}
		}
	}

	void IterateThroughNeighbors_cellList(size_t CenterParticle, double Rc, const std::function<void(const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)> callback, bool * pAbort = nullptr) const
	{
		this->IterateThroughNeighbors_cellList(this->GetRelativeCoordinates(CenterParticle), Rc, callback, pAbort);
	}



public:

	mutable unsigned long ParallelNumThread;
	//This class can work in one of the two modes:
	//Neighbor list mode: the class maintains a neighbor list, and update it if necessary using the cell list
	//Cell list mode (default): the class does not maintain a neighbor list
	//additionally, in cell list mode, there are two options:
	//UseSortedList=true : efficient when the basis vectors almost never change.
	//UseSortedList=false : efficient when the basis vectors change.
	PeriodicCellList()
	{
		//this->VicinityLatticesRc2=0;
		this->Dimension = 0;
		Index CellCount = 1;//need to calculate overall number of cells
		this->Cells.assign(CellCount, std::vector<size_t >());

		this->ApproximateCellSize = ::MaxDimension;
		this->SortedVicinityList = false;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		this->NeighborListMode = false;
		this->NeighborListMaxMovement2SinceBuilt = 0.0;
		this->NeighborListRange = 0.0;
		this->ParallelNumThread = 1;
		this->NeighborListCellDeformationCoeff = 10000;
		this->NeighborListLatticeShiftValid = false;
	}
	PeriodicCellList(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList = false)
	{
		//this->VicinityLatticesRc2=0;
		this->Dimension = Dimension;
		Index CellCount = 1;//need to calculate overall number of cells
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			GeometryVector & vec = BasisVectors[i];
			this->CellRank[i] = (signed long)(std::floor(std::sqrt(vec.Modulus2()) / CellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;

			//error check
			if (this->CellRank[i] <= 0)
			{
				std::cerr << "error in Configuration constructor: invalid cell rank. basis vector length squared=" << vec.Modulus2() << ", CellSize=" << CellSize << '\n';
				assert(false);
			}
			this->BasisVector[i] = vec;
			CellCount *= this->CellRank[i];
		}
		for (DimensionType i = this->Dimension; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = 0;
			this->BasisVector[i] = GeometryVector();
		}
		this->Cells.assign(CellCount, std::vector<size_t >());

		this->ApproximateCellSize = CellSize;
		this->SortedVicinityList = UseSortedList;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		this->NeighborListMode = false;
		this->NeighborListMaxMovement2SinceBuilt = 0.0;
		this->NeighborListRange = 0.0;
		this->ParallelNumThread = 1;
		this->NeighborListCellDeformationCoeff = 10000;
		this->NeighborListLatticeShiftValid = false;
	}
	PeriodicCellList(const PeriodicCellList<OtherCharacteristics> & source)
		: Dimension(source.Dimension),
		Cells(source.Cells), ApproximateCellSize(source.ApproximateCellSize),
		SortedVicinityList(source.SortedVicinityList), VicinityLatticesRc2(0.0),
		ParticleCartesians(source.ParticleCartesians), ParticleRelatives(source.ParticleRelatives),
		ParticleCharacteristics(source.ParticleCharacteristics)
	{
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->BasisVector[i] = source.BasisVector[i];
		}
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = source.CellRank[i];
		}
		this->ClearVicinityLatticeList();


		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		this->NeighborListMode = source.NeighborListMode;
		this->NeighborList = source.NeighborList;
		this->NeighborListMaxMovement2SinceBuilt = source.NeighborListMaxMovement2SinceBuilt;
		this->NeighborListStartingCartesians = source.NeighborListStartingCartesians;
		this->NeighborListRange = source.NeighborListRange;
		this->ParallelNumThread = source.ParallelNumThread;
		this->NeighborListCellDeformationCoeff = source.NeighborListCellDeformationCoeff;
		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->NeighborListStartingReciprocalBasis[i] = source.NeighborListStartingReciprocalBasis[i];
		this->NeighborListLatticeShiftValid = source.NeighborListLatticeShiftValid;
		this->NeighborListEnabledForThisParticle = source.NeighborListEnabledForThisParticle;
	}
	PeriodicCellList(std::istream & ifile, bool UseSortedList = false)
	{
		ifile.read((char *)(&this->Dimension), sizeof(this->Dimension));
		ifile.read((char *)(&ApproximateCellSize), sizeof(ApproximateCellSize));
		Index CellCount = 1;//need to calculate overall number of cells
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			this->BasisVector[i].ReadBinary(ifile, this->Dimension);
			GeometryVector & vec = this->BasisVector[i];
			this->CellRank[i] = (signed long)(std::floor(std::sqrt(vec.Modulus2()) / ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			//error check
			if (this->CellRank[i] <= 0)
			{
				std::cerr << "error in Configuration constructor: invalid cell rank. basis vector length squared=" << vec.Modulus2() << ", CellSize=" << ApproximateCellSize << '\n';
				assert(false);
			}
			CellCount *= this->CellRank[i];
		}
		for (DimensionType i = this->Dimension; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = 0;
			this->BasisVector[i] = GeometryVector();
		}
		this->Cells.assign(CellCount, std::vector<size_t >());

		this->ApproximateCellSize = ApproximateCellSize;
		this->SortedVicinityList = UseSortedList;
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		int nbrp;
		ifile.read((char *)(&nbrp), sizeof(nbrp));
		for (int i = 0; i<nbrp; i++)
		{
			GeometryVector temp;
			temp.ReadBinary(ifile, this->Dimension);
			OtherCharacteristics otemp;
			ifile.read((char *)(&otemp), sizeof(OtherCharacteristics));
			this->Insert(otemp, temp);
		}

		this->NeighborListMode = false;
		this->NeighborListMaxMovement2SinceBuilt = 0.0;
		this->NeighborListRange = 0.0;
		this->ParallelNumThread = 1;
	}
	template<typename add> PeriodicCellList(const PeriodicCellList<add> & source, const OtherCharacteristics & cha)
		: Dimension(source.Dimension),
		Cells(source.Cells), ApproximateCellSize(source.ApproximateCellSize),
		SortedVicinityList(source.SortedVicinityList), VicinityLatticesRc2(0.0),
		ParticleCartesians(source.ParticleCartesians), ParticleRelatives(source.ParticleRelatives),
		ParticleCharacteristics(source.NumParticle(), cha)
	{
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->BasisVector[i] = source.BasisVector[i];
		}
		for (DimensionType i = 0; i< ::MaxDimension; i++)
		{
			this->CellRank[i] = source.CellRank[i];
		}
		this->ClearVicinityLatticeList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		this->NeighborListMode = source.NeighborListMode;
		this->NeighborList = source.NeighborList;
		this->NeighborListMaxMovement2SinceBuilt = source.NeighborListMaxMovement2SinceBuilt;
		this->NeighborListStartingCartesians = source.NeighborListStartingCartesians;
		this->NeighborListRange = source.NeighborListRange;
		this->ParallelNumThread = source.ParallelNumThread;
		this->NeighborListCellDeformationCoeff = source.NeighborListCellDeformationCoeff;
		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->NeighborListStartingReciprocalBasis[i] = source.NeighborListStartingReciprocalBasis[i];
		this->NeighborListLatticeShiftValid = source.NeighborListLatticeShiftValid;
		this->NeighborListEnabledForThisParticle = source.NeighborListEnabledForThisParticle;
	}

	virtual ~PeriodicCellList()
	{
	}

	//switch to neighbor list mode
	void SwitchToNeighborListMode(const std::vector<char> & perParticleNeighborListEnabled = std::vector<char>()) const
	{
		if (!this->NeighborListMode)
		{
			this->NeighborListMode = true;
			this->NeighborListMaxMovement2SinceBuilt = 0.0;
			this->NeighborListRange = 0.0;
			this->NeighborListCellDeformationCoeff = 10000;
			this->NeighborListLatticeShiftValid = false;
			this->NeighborListStartingCartesians = this->ParticleCartesians;
			for (int i = 0; i < this->Dimension; i++)
				this->NeighborListStartingReciprocalBasis[i] = this->GetReciprocalBasisVector(i);

			if (perParticleNeighborListEnabled.size() == this->NumParticle())
				this->NeighborListEnabledForThisParticle = perParticleNeighborListEnabled;
			else
				this->NeighborListEnabledForThisParticle = std::vector<char>(this->NumParticle(), true);
		}
	}
	//swith off neighbor list
	void SwitchToCellListMode() const
	{
		if (this->NeighborListMode)
		{
			this->NeighborListMode = false;
			this->NeighborListMaxMovement2SinceBuilt = 0.0;
			this->NeighborListRange = 0.0;
			this->NeighborList.clear();
			this->NeighborListStartingCartesians.clear();
			this->RefineCellList();
			this->NeighborListEnabledForThisParticle.clear();
		}
	}
	bool UsingNeighborList() const
	{
		return this->NeighborListMode;
	}

	bool operator == (const PeriodicCellList<OtherCharacteristics> & right) const
	{
		if (right.Dimension != this->Dimension)
			return false;
		if (right.NumParticle() != this->NumParticle())
			return false;
		for (DimensionType d = 0; d<this->Dimension; d++)
			if (right.BasisVector[d] != this->BasisVector[d])
				return false;
		for (size_t i = 0; i<this->NumParticle(); i++)
			if (right.ParticleRelatives[i] != this->ParticleRelatives[i])
				return false;

		return true;
	}
	void ClearVicinityLatticeList(void) const
	{
		this->VicinityLattices.clear();
		this->VicinityLatticesRc2 = 0.0;
		this->VicinityLatticePartitions.clear();
		this->VicinityLatticePartitions.insert(std::make_pair((double)(0.0), (size_t)(0)));
	}
	void TryRefineBasisVectors_inner(int ilongest)
	{
		//find the other basis vector
		int iother;
		signed int movepart = 0;
		for (iother = 0; iother<this->Dimension; iother++)
		{
			if (iother == ilongest) continue;
			movepart = (signed int)(std::floor(0.5 + (this->BasisVector[ilongest].Dot(this->BasisVector[iother])) / (this->BasisVector[iother].Dot(this->BasisVector[iother]))));
			if (movepart != 0)
				break;
		}
		if (movepart == 0)
			return;

		this->BasisVector[ilongest] = this->BasisVector[ilongest] - static_cast<double>(movepart)*this->BasisVector[iother];
		for (size_t i = 0; i<this->ParticleRelatives.size(); i++)
		{
			this->ParticleRelatives[i].x[iother] += this->ParticleRelatives[i].x[ilongest] * movepart;
			this->PutToPeriodicBox(i);
			this->UpdateCartesianCoordinates(i);
		}
		this->ClearVicinityLatticeList();
		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;
		size_t NumCells = 1;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double length = std::sqrt(this->BasisVector[i].Modulus2());
			this->CellRank[i] = (signed long)(std::floor(length / this->ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			NumCells *= this->CellRank[i];
		}
		this->Cells.resize(NumCells);
		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();

		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->Cells[this->GetIndex(i)].push_back(i);
		}

		this->NeighborListRange = 0.0;
	}
	void TryRefineBasisVectors(void)
	{
		if (this->Dimension == 1)
			return;
		for (int i = 0; i<this->Dimension; i++)
			this->TryRefineBasisVectors_inner(i);
	}
	double PeriodicVolume(void) const
	{
		double result = ::Volume(this->BasisVector, this->Dimension);
		return result;
	}
	DimensionType GetDimension(void) const
	{
		return this->Dimension;
	}
	virtual void Insert(OtherCharacteristics cha, const GeometryVector & RelativeCoordinate)
	{
		//this->Particles.push_back(particle(this->Dimension, cha));
		this->ParticleCartesians.push_back(GeometryVector(this->Dimension));
		this->ParticleRelatives.push_back(GeometryVector(this->Dimension));
		this->ParticleCharacteristics.push_back(cha);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			assert(RelativeCoordinate.x[i] == RelativeCoordinate.x[i]);
			//this->Particles.back().RelativeCoordinates.x[i]=RelativeCoordinate.x[i]-std::floor(RelativeCoordinate.x[i]);
			//this->Particles.back().RelativeCoordinates.x[i]-=std::floor(this->Particles.back().RelativeCoordinates.x[i]);

			//NEED TO subtract the floor TWICE
			this->ParticleRelatives.back().x[i] = RelativeCoordinate.x[i] - std::floor(RelativeCoordinate.x[i]);
			this->ParticleRelatives.back().x[i] -= std::floor(this->ParticleRelatives.back().x[i]);
		}
		this->UpdateCartesianCoordinates(this->NumParticle() - 1);
		this->Cells[this->GetIndex(this->NumParticle() - 1)].push_back(this->NumParticle() - 1);

		this->NeighborListRange = 0.0;
	}
	void Insert(OtherCharacteristics cha, RandomGenerator & gen)
	{
		GeometryVector a(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			a.x[i] = gen.RandomDouble();
		this->Insert(cha, a);
	}
	void GetVicinityLattice(double Rc) const
	{
		size_t PrevVicinityLatticeNum = this->VicinityLattices.size();
		//prepare orthogonalized basis
		GeometryVector * ortho = new GeometryVector[this->Dimension];
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			ortho[i] = this->BasisVector[i];
			for (DimensionType j = 0; j<i; j++)
			{
				GeometryVector & e = ortho[j];
				ortho[i].MinusFrom(e*(ortho[i].Dot(e) / e.Dot(e)));
			}
		}

		this->_VicinityLatticeIteration(GeometryVector(this->Dimension), this->Dimension, GeometryVector(this->Dimension), ortho, Rc*Rc, this->VicinityLatticesRc2);

		delete[] ortho;

		if (this->SortedVicinityList)
			std::sort(this->VicinityLattices.begin() + PrevVicinityLatticeNum, this->VicinityLattices.end(), [](const VicinityLatticesMember & left, const VicinityLatticesMember & right) -> bool {return left.CartesianCoordinates.Modulus2()<right.CartesianCoordinates.Modulus2(); });

		this->VicinityLatticesRc2 = Rc*Rc;
	}
	void _PrintVicinityLattices(std::ostream & out) const
	{
		for (auto iter = this->VicinityLattices.begin(); iter<this->VicinityLattices.end(); iter++)
			out << (*iter).CartesianCoordinates << '\n';
	}
	size_t GetRandomParticle(RandomGenerator & gen, Index * SourceCell = nullptr) const//pick up a random particle, return it's pointer, if SourceCell is not null, also return the cell that the particle come from
	{
		size_t result;
		result = std::floor(gen.RandomDouble()*this->NumParticle());
		if (SourceCell != nullptr)
			*SourceCell = this->GetIndex(result);
		return result;
	}
	size_t NumParticle(void) const
	{
		return this->ParticleRelatives.size();
	}
	virtual OtherCharacteristics & GetCharacteristics(size_t Num)
	{
		return this->ParticleCharacteristics[Num];
	}
	virtual const OtherCharacteristics & GetCharacteristics(size_t Num) const
	{
		return this->ParticleCharacteristics[Num];
	}
	const GeometryVector & GetRelativeCoordinates(size_t Num) const
	{
		return this->ParticleRelatives[Num];
	}
	const GeometryVector & GetCartesianCoordinates(size_t Num) const
	{
		return this->ParticleCartesians[Num];
	}
	//return the length of the longest diagonal of the simulation box
	double MaxDiagonalLength(void) const
	{
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i]);
			else
				maxdist.MinusFrom(this->BasisVector[i]);
		}
		return std::sqrt(maxdist.Modulus2());
	}

	//allocate and generate internal cell lists, so that calling IterateThroughNeighbors with the same Rc is thread-safe
	void PrepareIterateThroughNeighbors_cellList(double Rc) const
	{
		//debug temp
		//std::cerr<<"IterateThroughNeighbors for Rc="<<Rc<<'\n';

		Rc = std::abs(Rc);
		double Rc2 = Rc*Rc;
		//find the true Rc, which consideres the distance from a particle to a lattice point
		GeometryVector maxdist(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			if (maxdist.Dot(this->BasisVector[i])>0)
				maxdist.AddFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
			else
				maxdist.MinusFrom(this->BasisVector[i] * (1.0 / this->CellRank[i]));
		}
		double newRc = Rc + std::sqrt(maxdist.Modulus2());
		double newRc2 = newRc*newRc;

		auto IterationEnd = this->VicinityLattices.end();
		//Get new VicinityLattice if necessary
		if (this->SortedVicinityList)
		{
			if ((newRc2)>this->VicinityLatticesRc2)
				this->GetVicinityLattice(newRc*1.1);
			IterationEnd = this->VicinityLattices.end();
		}
		else
		{
			auto ub = this->VicinityLatticePartitions.lower_bound(newRc2);
			auto lb = ub;
			if (this->VicinityLatticePartitions.size()>1)
				lb--;
			if ((newRc2)>this->VicinityLatticesRc2)
			{
				this->GetVicinityLattice(newRc);
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, this->VicinityLattices.size()));
				IterationEnd = this->VicinityLattices.end();
			}
			else if (ub->first == newRc2)
			{
				IterationEnd = ub->second + this->VicinityLattices.begin();
			}
			else if (ub->second - lb->second>10)
			{
				auto titer = std::partition(this->VicinityLattices.begin() + lb->second, this->VicinityLattices.begin() + ub->second,
					[&newRc2](const VicinityLatticesMember & a) ->bool
				{
					return a.CartesianCoordinates.Modulus2()<newRc2;
				});
				this->VicinityLatticePartitions.insert(std::make_pair(newRc2, titer - this->VicinityLattices.begin()));
				IterationEnd = titer;
			}
			else
				IterationEnd = ub->second + this->VicinityLattices.begin();
		}
	}

	//iterate through all neighbors centered at RelativeCoordinate with radius Rc.
	//parameters of callback function:
	//shift : vector from starting point to neighbor. LatticeShift : the part of shift that is because the starting point and neighbor is not in the same periodic cell
	//PeriodicShift : same as LatticeShift, expressed in multiples of each lattice vectors. Sourceparticle : particle pointer to neighbor.
	//pAbort : if it is not nullptr and points to true, the iteration will be stopped immediately
	//this version cannot use neighbor list, since the iteration center is not necessarily a particle
	void IterateThroughNeighbors(const GeometryVector & RelativeCoordinate, double Rc, const std::function<void(const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)> callback, bool * pAbort = nullptr) const
	{
		if (this->NeighborListMode)
		{
			std::cerr << "Warning in PeriodicCellList::IterateThroughNeighbors : Iterating through neighbors of an arbitrary center cannot be performed in neighbor list mode!\n";
			std::cerr << "Switching to cell list mode!\n";
			this->SwitchToCellListMode();
		}

		this->IterateThroughNeighbors_cellList(RelativeCoordinate, Rc, callback, pAbort);
	}

	//is the neighbor list ready for an iteration of range Rc?
	bool NeighborListReady(double Rc) const
	{
		return (this->NeighborListCellDeformationCoeff*Rc + 2 * std::sqrt(NeighborListMaxMovement2SinceBuilt)) < NeighborListRange;
	}
	void BuildNeighborList(double Range) const
	{
		if (this->NeighborListMode == false)
		{
			std::cerr << "Error in BuildNeighborList : not in neighbor list mode!\n";
			assert(false);
		}
		this->NeighborList = std::vector< std::vector<NeighborListElementType> >(this->NumParticle(), std::vector<NeighborListElementType>());
		this->RefineCellList();
		this->PrepareIterateThroughNeighbors_cellList(Range);
#pragma omp parallel for num_threads(this->ParallelNumThread)
		for (int i = 0; i < this->NumParticle(); i++)
		{
			if (NeighborListEnabledForThisParticle[i])
			{
				std::vector<NeighborListElementType> & l = this->NeighborList[i];
				auto func = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
				{
					NeighborListElementType t;
					t.SourceParticle = Sourceparticle;
					std::memcpy(t.PeriodicShift, PeriodicShift, sizeof(const signed long)*::MaxDimension);
					t.LatticeShift = LatticeShift;
					l.push_back(t);
				};
				this->IterateThroughNeighbors_cellList(i, Range, func);
			}
		}
		this->NeighborListStartingCartesians = this->ParticleCartesians;
		this->NeighborListMaxMovement2SinceBuilt = 0.0;
		this->NeighborListRange = Range;

		for (DimensionType i = 0; i < this->Dimension; i++)
			this->NeighborListStartingReciprocalBasis[i] = this->GetReciprocalBasisVector(i);
		this->NeighborListCellDeformationCoeff = 1.0;

		this->NeighborListLatticeShiftValid = true;

		//debug temp
		//std::cout << "neighbor list built\n";
	}
	void PrepareIterateThroughNeighbors(double Rc) const
	{
		if (this->NeighborListMode)
		{
			if (!NeighborListReady(Rc))
				BuildNeighborList(1.1*Rc);
		}
		else
			PrepareIterateThroughNeighbors_cellList(Rc);
	}
	void IterateThroughNeighbors(size_t CenterParticle, double Rc, const std::function<void(const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)> callback, bool * pAbort = nullptr) const
	{
		double Rc2 = Rc*Rc;
		if (this->NeighborListMode && this->NeighborListEnabledForThisParticle[CenterParticle])
		{
			if (!NeighborListReady(Rc))
				BuildNeighborList(1.1*Rc);

			if (!this->NeighborListLatticeShiftValid)
			{
				//re-calculate lattice shift
				for (auto iter = this->NeighborList.begin(); iter != this->NeighborList.end(); ++iter)
				{
					std::vector<NeighborListElementType> & l = *iter;
					for (auto iter = l.begin(); iter != l.end(); ++iter)
					{
						GeometryVector latticeShift(this->Dimension);
						for (DimensionType d = 0; d < this->Dimension; d++)
							latticeShift.AddFrom(iter->PeriodicShift[d] * this->BasisVector[d]);
						iter->LatticeShift = latticeShift;
					}
				}
				this->NeighborListLatticeShiftValid = true;
			}

			//do iteration
			std::vector<NeighborListElementType> & l = this->NeighborList[CenterParticle];
			GeometryVector car = this->GetCartesianCoordinates(CenterParticle);
			for (auto iter = l.begin(); iter != l.end(); ++iter)
			{
				GeometryVector shift = this->GetCartesianCoordinates(iter->SourceParticle) - car + iter->LatticeShift;
				if (shift.Modulus2() < Rc2)
					callback(shift, iter->LatticeShift, iter->PeriodicShift, iter->SourceParticle);

				if (pAbort != nullptr)
					if (*pAbort == true)
						return;
			}
		}
		else if (this->NeighborListMode)
			std::cerr << "Error in IterateThroughNeighbors : neighbor list enabled, but disabled for the particular center particle\n";
		else
			IterateThroughNeighbors_cellList(CenterParticle, Rc, callback, pAbort);
	}

	void MoveParticle(size_t Num, const GeometryVector & newRelativeCoordinate)
	{
		if (this->NeighborListMode)
		{
			this->ParticleRelatives[Num] = newRelativeCoordinate;
			for (DimensionType i = 0; i < this->Dimension; i++)
			{
				this->ParticleRelatives[Num].x[i] -= std::floor(this->ParticleRelatives[Num].x[i]);
				this->ParticleRelatives[Num].x[i] -= std::floor(this->ParticleRelatives[Num].x[i]);
			}
			this->UpdateCartesianCoordinates(Num);
			double movement2 = (this->ParticleCartesians[Num] - this->NeighborListStartingCartesians[Num]).Modulus2();
			this->NeighborListMaxMovement2SinceBuilt = std::max(this->NeighborListMaxMovement2SinceBuilt, movement2);
		}
		else
		{
			//particle * a=&this->Particles[Num];
			Index fromIndex;
			fromIndex = this->GetIndex(Num);
			std::vector<size_t > & currentCell = this->Cells[fromIndex];
			std::vector<size_t >::iterator itera = currentCell.begin();
			for (;; itera++)
			{
				assert(itera != currentCell.end());
				if (*itera == Num)
					break;
			}

			this->ParticleRelatives[Num] = newRelativeCoordinate;
			for (DimensionType i = 0; i < this->Dimension; i++)
			{
				this->ParticleRelatives[Num].x[i] -= std::floor(this->ParticleRelatives[Num].x[i]);
				this->ParticleRelatives[Num].x[i] -= std::floor(this->ParticleRelatives[Num].x[i]);
			}
			this->UpdateCartesianCoordinates(Num);
			this->UpdateCellList(itera, currentCell);
		}
	}
	void DeleteParticle(size_t Num)
	{
		//particle * a=&this->Particles[Num];
		Index afromIndex = this->GetIndex(Num);
		std::vector<size_t > & acurrentCell = this->Cells[afromIndex];
		std::vector<size_t >::iterator itera = acurrentCell.begin();
		for (;; itera++)
		{
			assert(itera != acurrentCell.end());
			if (*itera == Num)
				break;
		}
		if (Num == this->NumParticle() - 1)
		{
			this->ParticleCartesians.pop_back();
			this->ParticleRelatives.pop_back();
			this->ParticleCharacteristics.pop_back();
			if (itera == acurrentCell.end() - 1)
				acurrentCell.pop_back();
			else
			{
				(*itera) = acurrentCell.back();
				acurrentCell.pop_back();
			}
			return;
		}

		//particle * back=&this->Particles.back();
		size_t backi = this->NumParticle() - 1;
		Index bfromIndex = this->GetIndex(backi);
		std::vector<size_t > & bcurrentCell = this->Cells[bfromIndex];
		std::vector<size_t >::iterator iterb = bcurrentCell.begin();
		for (;; iterb++)
		{
			assert(iterb != bcurrentCell.end());
			if (*iterb == (backi))
				break;
		}

		//(*a)=this->Particles.back();
		//this->Particles.pop_back();
		this->ParticleCartesians[Num] = this->ParticleCartesians.back();
		this->ParticleRelatives[Num] = this->ParticleRelatives.back();
		this->ParticleCharacteristics[Num] = this->ParticleCharacteristics.back();
		this->ParticleCartesians.pop_back();
		this->ParticleRelatives.pop_back();
		this->ParticleCharacteristics.pop_back();

		(*iterb) = Num;
		if (itera == acurrentCell.end() - 1)
			acurrentCell.pop_back();
		else
		{
			(*itera) = acurrentCell.back();
			acurrentCell.pop_back();
		}

		this->NeighborListRange = 0.0;
	}
	GeometryVector GetBasisVector(::DimensionType i) const
	{
		return this->BasisVector[i];
	}
	GeometryVector GetReciprocalBasisVector(::DimensionType i) const
	{
		if (this->ReciprocalBasisVectorValid[i])
			return this->ReciprocalBasisVector[i];
		else
		{

			GeometryVector result = this->BasisVector[i];

			/////////////////////////////////////////////
			std::vector<GeometryVector> ortho;
			ortho.reserve(this->Dimension);
			for (DimensionType j = 0; j<this->Dimension; j++)
			{
				if (j == i) continue;
				GeometryVector temp = this->BasisVector[j];
				for (auto iter = ortho.begin(); iter != ortho.end(); iter++)
					temp.MinusFrom((*iter)*((temp.Dot(*iter)) / (iter->Modulus2())));
				ortho.push_back(temp);
			}
			////////////////////////////////////////////

			for (auto iter = ortho.begin(); iter != ortho.end(); iter++)
				result.MinusFrom((*iter)*((iter->Dot(result)) / (iter->Modulus2())));

			double product = result.Dot(this->BasisVector[i]);
			result.MultiplyFrom(2 * ::pi / product);
			this->ReciprocalBasisVector[i] = result;
			this->ReciprocalBasisVectorValid[i] = true;
			return result;
		}
	}
	void RefineCellList(void) const
	{
		this->ClearVicinityLatticeList();
		size_t NumCells = 1;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double length = std::sqrt(this->BasisVector[i].Modulus2());
			this->CellRank[i] = (signed long)(std::floor(length / this->ApproximateCellSize));
			if (this->CellRank[i] == 0)
				this->CellRank[i] = 1;
			NumCells *= this->CellRank[i];
		}
		this->Cells.resize(NumCells);
		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();

		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->Cells[this->GetIndex(i)].push_back(i);
		}
	}
	double GetCellSize(void)
	{
		return this->ApproximateCellSize;
	}
	void SetSortedCellList(bool SortedCellList)
	{
		if (SortedCellList != this->SortedVicinityList)
			this->ClearVicinityLatticeList();
		this->SortedVicinityList = SortedCellList;
	}
	void SetCellSize(double CellSize)
	{
		this->ApproximateCellSize = CellSize;
		bool NeedRefine = false;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double celllength = std::sqrt(this->BasisVector[i].Modulus2()) / this->CellRank[i];
			if (celllength > 2 * this->ApproximateCellSize)
				NeedRefine = true;
			else if (celllength < 0.5 * this->ApproximateCellSize && this->CellRank[i] > 1)
				NeedRefine = true;
		}
		if (NeedRefine)
			this->RefineCellList();
	}
	void ChangeBasisVector(const GeometryVector NewBasisVectors[])
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->BasisVector[i] = NewBasisVectors[i];
		this->ClearVicinityLatticeList();
		for (size_t i = 0; i<this->NumParticle(); i++)
		{
			this->UpdateCartesianCoordinates(i);
		}

		bool NeedRefine = false;
		for (DimensionType i = 0; i<this->Dimension; i++)
		{
			double celllength = std::sqrt(this->BasisVector[i].Modulus2()) / this->CellRank[i];
			if (celllength > 2 * this->ApproximateCellSize)
				NeedRefine = true;
			else if (celllength < 0.5 * this->ApproximateCellSize && this->CellRank[i] > 1)
				NeedRefine = true;
		}
		if (NeedRefine)
			this->RefineCellList();

		for (DimensionType i = 0; i< ::MaxDimension; i++)
			this->ReciprocalBasisVectorValid[i] = false;

		if (this->NeighborListMode)
			this->NeighborListCellDeformationCoeff = this->getLambdaMaxOfDeformation();

		this->NeighborListLatticeShiftValid = false;
	}
	void RemoveParticles(void)//remove all Particles
	{
		//this->Particles.clear();
		this->ParticleCartesians.clear();
		this->ParticleCharacteristics.clear();
		this->ParticleRelatives.clear();

		for (auto iter = this->Cells.begin(); iter != this->Cells.end(); iter++)
			iter->clear();

		this->NeighborListRange = 0.0;
	}
	void Rescale(double Factor)
	{
		this->ApproximateCellSize *= Factor;
		std::vector<GeometryVector> basis;
		for (DimensionType i = 0; i<this->Dimension; i++)
			basis.push_back(this->BasisVector[i] * Factor);
		this->ChangeBasisVector(&basis[0]);
	}
	void Resize(double NewVolume)
	{
		double factor = std::pow(NewVolume / this->PeriodicVolume(), 1.0 / this->Dimension);
		this->Rescale(factor);
	}
	GeometryVector CartesianCoord2RelativeCoord(const GeometryVector & CartesianCoord) const
	{
		GeometryVector result(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			result.x[i] = CartesianCoord.Dot(this->GetReciprocalBasisVector(i));

		return result*(0.5 / pi);
	}
	GeometryVector RelativeCoord2CartesianCoord(const GeometryVector & RelativeCoord) const
	{
		GeometryVector result(this->Dimension);
		for (DimensionType i = 0; i<this->Dimension; i++)
			result.AddFrom((this->BasisVector[i])*(RelativeCoord.x[i]));

		return result;
	}
	void RelativeCoordToMinimumImage(GeometryVector & a) const
	{
		for (DimensionType i = 0; i<this->Dimension; i++)
			a.x[i] -= std::floor(a.x[i]+0.5);
	}

	//write the configuration in binary, this saves space compared to Output()
	// there is no such thing as void ReadBinary(std::istream & ifile)! Use the constructor instead!
	void WriteBinary(std::ostream & ofile) const
	{
		ofile.write((char *)(&this->Dimension), sizeof(this->Dimension));
		ofile.write((char *)(&ApproximateCellSize), sizeof(ApproximateCellSize));
		for (DimensionType i = 0; i<this->Dimension; i++)
			this->BasisVector[i].WriteBinary(ofile, this->Dimension);
		int nbrp = this->NumParticle();
		ofile.write((char *)(&nbrp), sizeof(nbrp));
		for (int i = 0; i<nbrp; i++)
		{
			this->ParticleRelatives[i].WriteBinary(ofile, this->Dimension);
			ofile.write((char *)(&this->ParticleCharacteristics[i]), sizeof(OtherCharacteristics));
		}
	}

	//rotate the current structure so that
	//the 1st basis vector is (c11, 0, 0, ...)
	//the 2nd basis vector is (c21, c22, 0, ...)
	// ......
	//AND the basis vector length l1>=l2>=l3...
	void StandardizeOrientation(void)
	{
		::StandardizeOrientation(this->BasisVector, this->Dimension);
		this->ChangeBasisVector(this->BasisVector);
	}
	//PeriodicCellList< of anything > can be converted to PeriocicCellList<Empty>
	operator PeriodicCellList<Empty> const () const
	{
		return PeriodicCellList<Empty>(*this, Empty());
	}

	//return the distance from fromRelative to the closest particle
	GeometryVector NearestParticleVector(const GeometryVector & fromRelative) const
	{
		double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension())*1.0;
		std::vector<GeometryVector> neighbors;
		double l = TypicalLength;
		while (neighbors.size() < 1)
		{
			neighbors.clear();
			(*this).IterateThroughNeighbors(fromRelative, l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l *= 1.5;
		}

		std::partial_sort(neighbors.begin(), neighbors.begin() + 1, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2() < right.Modulus2(); });
		return neighbors[0];
	}
	//return the distance from fromRelative to the closest particle
	double NearestParticleDistance(const GeometryVector & fromRelative) const
	{
		return std::sqrt(NearestParticleVector(fromRelative).Modulus2());
	}

	//return the distance from particle No. i to the closest particle
	GeometryVector NearestParticleVector(size_t i) const
	{
		double TypicalLength = std::pow((*this).PeriodicVolume() / (*this).NumParticle(), 1.0 / (*this).GetDimension())*1.0;
		std::vector<GeometryVector> neighbors;
		double l = TypicalLength;
		while (neighbors.size() < 2)
		{
			neighbors.clear();
			(*this).IterateThroughNeighbors(i, l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l *= 1.5;
		}

		std::partial_sort(neighbors.begin(), neighbors.begin() + 2, neighbors.end(), [](const GeometryVector & left, const GeometryVector & right)->bool{return left.Modulus2() < right.Modulus2(); });
		return neighbors[1].Modulus2();
	}
	double NearestParticleDistance(size_t i) const
	{
		return std::sqrt(NearestParticleVector(i).Modulus2());
	}
};


//overlap all the configurations in that std::vector
//require that all the configurations in the std::vector have the same basis vectors
template<typename a>PeriodicCellList<a> CombinedList(const std::vector<PeriodicCellList<a> > & components)
{
	if(components.size()==0)
		return PeriodicCellList<a>();
	PeriodicCellList<a> result = components[0];
	for(size_t i=1; i<components.size(); i++)
	{
		const PeriodicCellList<a> & list = components[i];
		for(size_t j=0; j<list.NumParticle(); j++)
			result.Insert(list.GetCharacteristics(j), list.GetRelativeCoordinates(j));
	}
	return result;
}
//config is a CombinedList of m PeriodicCellLists, extract the nth list
template<typename a>PeriodicCellList<a> SeparateList(const PeriodicCellList<a> & config, size_t m, size_t n)
{
	size_t nt=config.NumParticle();
	if(nt%m!=0)
	{
		std::cerr<<"Error in SeparateList : the total number of particles is not multiple of the number of configurations!\n";
		assert(false);
	}
	size_t np=nt/m;
	size_t start=n*np;
	size_t end=start+np;
	PeriodicCellList<a> result = config;
	result.RemoveParticles();
	for(size_t i=start; i<end; i++)
	{
		result.Insert(config.GetCharacteristics(i), config.GetRelativeCoordinates(i));
	}
	return result;
}

struct AtomInfo
{
	char name[4];
	double radius;
	AtomInfo()
	{
		this->name[0]='\0';
		this->name[3]='\0';
		this->radius = 0.0;
	}
	AtomInfo(const char * n)
	{
		std::strncpy(this->name, n, 3);
		this->name[3] = '\0';
		this->radius = 0.0;
	}
	AtomInfo(const char * n, double r)
	{
		std::strncpy(this->name, n, 3);
		this->name[3] = '\0';
		this->radius = r;
	}
	bool operator == (const AtomInfo & right) const
	{
		return std::strncmp(name, right.name, 3) == 0;
	}
	bool operator != (const AtomInfo & right) const
	{
		return std::strncmp(name, right.name, 3) != 0;
	}
	bool operator < (const AtomInfo & right) const
	{
		return std::strncmp(name, right.name, 3) < 0;
	}
	operator const char*() const
	{
		return name;
	}
};



typedef PeriodicCellList<AtomInfo> Configuration;

void Output(std::string prefix, const Configuration & List);
void Output(std::ostream & out, const Configuration & List);
inline void WriteBinary(std::ostream & out, const Configuration & List)
{
	List.WriteBinary(out);
}
inline void WriteBinary(std::string prefix, const Configuration & List)
{
	prefix+=std::string(".Config");
	std::fstream ofile(prefix.c_str(), std::fstream::out | std::fstream::binary);
	WriteBinary(ofile, List);
}
inline Configuration ReadBinary(std::istream & ifile)
{
	return Configuration(ifile);
}
inline Configuration ReadBinary(std::string prefix)
{
	prefix+=std::string(".Config");
	std::fstream ifile(prefix.c_str(), std::fstream::in | std::fstream::binary);
	return ReadBinary(ifile);
}

Configuration ReadPos(std::string prefix);
Configuration ReadPos(std::istream & ifile);
Configuration ReadHoomdXml(std::istream & ifile);
Configuration ReadHoomdXml(std::string prefix);
Configuration ReadStealthOutput(std::istream & file, double SideLength);
Configuration ReadStealthOutput(std::istream & file);
Configuration MultiplicateStructure(const Configuration & src, size_t NumberInEachSide);
Configuration MultiplicateStructure(const Configuration & src, const std::vector<size_t> & NumberInEachSide);
Configuration ReadGro(std::istream & file);
Configuration GetUnitCubicBox(DimensionType d, double CellSize=0.1);

//generate Poisson configuration of number density 1, in a cubic box
Configuration GeneratePoisson(DimensionType d, size_t N, RandomGenerator & gen);
Configuration ReadCoordinate(std::istream & file);
Configuration ReadCoordinate(std::istream & file, double SideLength);
Configuration ReadCartesianCoordinate(std::istream & file);
Configuration ReadCartesianCoordinate(std::istream & file, double SideLength, DimensionType dimension=0);

// a pack of multiple configurations, recorded in files Prefix.ConfigPack and Prefix.ConfigPackIndex
// .ConfigPack file contains configurations in binary format, one by one
// .ConfigPackIndex has the format:
// long long Number of Config, long long location of config0, long long location of config1, ...
class ConfigurationPack
{
private:
	std::string IndexName, PackName;
	long long NConfig;
public:
	ConfigurationPack()
	{
		this->NConfig=0;
	}
	ConfigurationPack(const std::string & prefix)
	{
		this->Open(prefix);
	}
	void Open(const std::string & prefix);
	long long NumConfig(void) const;
	void AddConfig(const Configuration & c);
	Configuration GetConfig(long long i) const;

	//discard all configuurations that are already in the file
	void Clear(void);
};

#include "DisplaySpheres.h"
//if the two functions below is not customizable enough, 
//then use this function to convert a configuration to vectors of spheres and lines,
//do the desired customization
//and then use DisplaySpheres/DisplaySpheres_Movie
void ConfigurationToSpheres(const Configuration & list, std::vector<sphere> & spheres, std::vector<line> & lines, double Radius = 0.0, double BondLength = 0.0, const std::vector<GeometryVector> & displacement = std::vector<GeometryVector>());

//use these functions to display configurations in a new window.
void DisplayConfiguration(const Configuration & a, double BondLength=0.0, double Radius=0.0);

void Display3DConfigurationMovie(ConfigurationPack & pk, double BondLength = 0.0, double Radius = 0.0);

void Display3DStructureFactor(const Configuration& a, double KCutOff, double KSphereRadius=0.0);

//polydisperse sphere packing
class SpherePacking : public PeriodicCellList<double>
{
protected:
	mutable double MaxRadius;//=0 indicates this variable is not valid. call UpdateMaxRadius to update it
public:
	SpherePacking() : PeriodicCellList(), MaxRadius(0.0)
	{
	}
	SpherePacking(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList=false);
	template<typename add> SpherePacking(const PeriodicCellList<add> & source, const double & radius) : PeriodicCellList(source, radius)
	{
		this->MaxRadius=radius;
	}
	virtual double & GetCharacteristics(size_t Num)
	{
		this->MaxRadius = 0.0;
		return this->ParticleCharacteristics[Num];
	}
	virtual const double & GetCharacteristics(size_t Num) const
	{
		return this->ParticleCharacteristics[Num];
	}

	SpherePacking(std::istream & ifile, bool UseSortedList=false) : PeriodicCellList<double>(ifile, UseSortedList)
	{
		this->MaxRadius=0.0;
		for(int i=0; i<this->NumParticle(); i++)
			this->MaxRadius=std::max(this->MaxRadius, this->GetCharacteristics(i));
	}
	void UpdateMaxRadius(void) const
	{
		if (this->MaxRadius == 0.0)
		{
			for (int i = 0; i < this->NumParticle(); i++)
				this->MaxRadius = std::max(this->MaxRadius, this->ParticleCharacteristics[i]);
		}
	}
	double GetMaxRadius(void) const
	{
		UpdateMaxRadius();
		return this->MaxRadius;
	}
	virtual void Insert(double radius, const GeometryVector & RelativeCoordinate);

	//check if the voxel is COMPLETELY occupied
	//no need to distinguish partial availability and full availability
	//relative coord and halfsize 
	//Radius : cartesian radius of the upcoming sphere
	bool CheckVoxelOverlap(const GeometryVector & cr, double halfsize, double Radius) const;

	//returns 0 if voxel is completely empty
	//2 if voxel is fully occupied
	//1 otherwise
	//relative coord and halfsize 
	//Radius : cartesian radius of the upcoming sphere
	int CheckVoxelOverlap2(const GeometryVector & cr, double halfsize, double Radius) const;


	//relative coord
	//Radius : cartesian radius of the upcoming sphere
	bool CheckOverlap(const GeometryVector & cr, double Radius) const;

	//check if there are overlapping spheres
	bool CheckOverlap(void) const
	{
		UpdateMaxRadius();
		for (size_t i = 0; i < this->NumParticle(); i++)
		{
			bool Overlap = false;
			const std::vector<double> * pSphereRadii = &this->ParticleCharacteristics;
			double Radius = this->GetCharacteristics(i);
			this->IterateThroughNeighbors(i, this->GetCharacteristics(i) + MaxRadius, [&Overlap, &Radius, &pSphereRadii, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
			{
				if (Sourceparticle != i || LatticeShift.Modulus2() != 0)
				{
					double OverlapRadius = Radius + pSphereRadii->at(Sourceparticle);
					if (shift.Modulus2() < OverlapRadius*OverlapRadius)
						Overlap = true;
				}
			}, &Overlap);
			if (Overlap)
				return true;
		}
		return false;
	}

	double PackingFraction(void) const
	{
		double v1=0.0;
		double d=this->GetDimension();
		for(int i=0; i<this->NumParticle(); i++)
			v1+= ::HyperSphere_Volume(d, this->GetCharacteristics(i));
		return v1/this->PeriodicVolume();
	}
};


SpherePacking ReadTJOutput(std::istream & inFile);

SpherePacking ReadTJOutput(const std::string & prefix);


void WriteTJOutput(std::fstream & ofile, const SpherePacking & pk);
void WriteTJOutput(const std::string & prefix, const SpherePacking & pk);

//VoxelType should be either char or bool. Char is faster, while bool can work with smaller memory.
typedef char VoxelType;
//pixelizing the configuration to MeshSide pixes in each side
void DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> & result);
inline std::vector<VoxelType> DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide)
{
	std::vector<VoxelType> result;
	DigitizeConfiguration(c, radius, MeshSide, result);
	return result;
}
//calculate volume fraction if each particle is replaced with a hypersphere of certain radius
//do this by pixelizing the configuration to MeshSide pixes in each side
double Volume(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> * pbuffer = nullptr);


//return the minimum D such that if each particle is replaced with a sphere with diameter D, particles a and b will be connected.
double MinConnectDiameter(const Configuration & c, size_t a, size_t b);
//return the minimum D such that if each particle is replaced with a sphere with diameter D, 
//particle n will connect to its periodic image in direction dir.
double PercolationDiameter(const Configuration & c, size_t n, DimensionType dir);


//Replace each particle in c with a hypersphere with radius,
//what's the cluster size InitParticle is in?
//FoundPeriodicImages: whether or not found more than 1 periodic images of a particle, 
//which means the cluster size is actually infinite rather than the returned value.
size_t ClusterSize(const Configuration & c, size_t InitParticle, double radius, bool & FoundPeriodicImages);



#include <complex>
//return Psi_6 of particle j in Config
//in 2D, Psi_6=1/N_{neighbors}( \\Sum_{neighbors} exp(6i \\theta) ).
//in 3D, an additional parameter, m, has to be provided. Psi_6=1/N_{neighbors}( \\Sum_{neighbors} exp(m * i \\phi)*P_6^m(cos(\\theta)) ), where P_l^m is the normalized associated Legendre polynomial.
//definition of neighbors depend on NeighborType
//0: six nearest neighbors in 2D
//1: neighbor in Voronoi diagram
//where theta is the angle between the bond and some reference direction
std::complex<double> Psi6(const PeriodicCellList<Empty> & Config, size_t j, int NeighborType = 0, signed int m=0);


//This class initializes by assigining each particle in a PeriodicCellList to its own cluster
//Then, as the user call member function Proceed(double r), particles within distance r are merged to a single cluster
//The user can then query various member objects of this class for various aspects of clusters:
//member clusters is a list of all clusters. Each cluster is represented as a list of its members' index
//affiliation[i] is the cluster particle i currently belong to
class ClusterMerger
{
public:
	typedef std::list<size_t> cluster;
	struct Bond
	{
		double l;
		size_t i, j;
	};
private:
	//forbid copy
	ClusterMerger(const ClusterMerger & src);
	std::vector<Bond> vBond;
	std::vector<Bond>::iterator next;
	size_t privateN;
	std::list<cluster> _clusters;
	std::vector<std::list<cluster>::iterator> _affiliation;
	size_t _LargestClusterSize;
public:
	const size_t & N=privateN;

//clusters is a list of all clusters. Each cluster is represented as a list of its members' index
	const std::list<cluster> & clusters=_clusters;

//affiliation[i] is the cluster particle i currently belong to
	const std::vector<std::list<cluster>::iterator> & affiliation=_affiliation;
	const size_t & LargestClusterSize=_LargestClusterSize;
	//currently, all bonds within distance currentL has been merged
	double currentL;


public:
	//manually specify distances between particles
	ClusterMerger(size_t N, const std::vector<Bond> & Bonds) : privateN(N), vBond(Bonds)
	{
		for (size_t i = 0; i < privateN; i++)
		{
			//at the beginning, each particle consist its own cluster
			cluster temp;
			temp.insert(temp.begin(), i);
			_clusters.insert(_clusters.begin(), temp);
			_affiliation.push_back(_clusters.begin());
		}
		std::sort(vBond.begin(), vBond.end(), [](const Bond & left, const Bond & right) ->bool
		{
			return left.l < right.l;
		});
		_LargestClusterSize = 1;
		next = vBond.begin();
	}

	//automatically find distances between particles
	ClusterMerger(const PeriodicCellList<Empty> & c, double Rmax)
	{
		privateN = c.NumParticle();
		for (size_t i = 0; i < privateN; i++)
		{
			//at the beginning, each particle consist its own cluster
			cluster temp;
			temp.insert(temp.begin(), i);
			_clusters.insert(_clusters.begin(), temp);
			_affiliation.push_back(_clusters.begin());
			std::vector<Bond> & vb = vBond;
			//generate a vector of bonds for future computation
			c.IterateThroughNeighbors(i, Rmax, [&vb, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				if (i < Sourceparticle)
				{
					Bond temp;
					temp.i = i;
					temp.j = Sourceparticle;
					temp.l = std::sqrt(shift.Modulus2());
					vb.push_back(temp);
				}
			});
		}
		if (vBond.size() == 0)
		{
			std::cerr << "Error in ClusterMerger constructor : no pair distance found!\n";
			return;
		}
		std::sort(vBond.begin(), vBond.end(), [](const Bond & left, const Bond & right) ->bool
		{
			return left.l < right.l;
		});
		_LargestClusterSize = 1;
		next = vBond.begin();
		currentL = next->l;
	}
	void Proceed(std::function<void(size_t, size_t)> mergeCallBack = [](size_t i, size_t j){})
	{
		std::list<cluster>::iterator pi = _affiliation[next->i], pj = _affiliation[next->j];
		if (pi != pj)
		{
			//merge the two _clusters
			cluster & cj = *_affiliation[next->j];
			cluster & ci = *_affiliation[next->i];
			//1. change the _affiliation of all particles in the same cluster as j
			for (auto iter = cj.begin(); iter != cj.end(); ++iter)
			{
				_affiliation[*iter] = _affiliation[next->i];
				mergeCallBack(*iter, next->i);
			}
			//2. merge the two _clusters
			ci.splice(ci.begin(), cj);
			_clusters.erase(pj);
			_LargestClusterSize = std::max(_LargestClusterSize, ci.size());
		}
		currentL = next->l;
		next++;
	}
	//merge _clusters until distances of any cluster is larger than r
	void Proceed(double r, std::function<void(size_t, size_t)> mergeCallBack = [](size_t i, size_t j){})
	{
		while (next != vBond.end() && next->l < r)
		{
			this->Proceed(mergeCallBack);
		}
	}
};


//display g1(r) of a 3D ensemble
//g1(r) is displayed using spheres in a 3D grid. Sphere volume is proportional to g1.
//sphere color correspond to ln(g1+1), when this quantity increases from 0 to its maximum,
//color changes from black to blue, cyan, yellow, and red
void DisplayG1_3D(ConfigurationPack pk, size_t NGridPerSide);


//test if two configurations are equivalent
//works by identifying the closest pair, second closest pair, ..., dth closest pair,
//and use these pairs to find rotation/inversion matrix and the translation vector.
//Therefore, this function works well only if the configurations are disordered.
bool ConfigurationEquivalent(const Configuration & a, const Configuration & b);


#endif
