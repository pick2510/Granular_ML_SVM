#include "RandomSequentialAddition.h"

struct SpherePacking_TempStruct
{
	signed char c[MaxDimension];
	float t;
};
const double SystemSize=1;
//convert double * to GeometryVector
GeometryVector VoxelList::Coord2Vect(const double * coord)
{
	GeometryVector cr(this->dimension);
	std::memcpy(cr.x, coord, this->dimension*sizeof(double));
	return cr;
}

void VoxelList::getsplitvoxellist(void)
{
	int max=(int)(std::pow((double)(2),this->dimension)+0.1);
	this->nbrsplitvoxel=max;
	this->splitvoxellist= new signed char [max*MaxDimension];
	signed char increment=2;
	signed char finish=1;
	signed char coord[MaxDimension];
	int i;
	std::vector<SpherePacking_TempStruct> result;
	result.resize(max);
	std::vector<SpherePacking_TempStruct>::iterator iter=result.begin();
	//initialize
	for(i=0; i<this->dimension; i++)
		coord[i]=-1;
	//main loop body
start:
	for(i=0; i<this->dimension; i++)
	{
		iter->c[i]=coord[i];
	}
	iter++;
	//loop end
	coord[0]+=increment;
	for(i=0; i<this->dimension-1; i++)
		if(coord[i]>finish)
		{
			coord[i]=-1;
			coord[i+1]+=increment;
		}
		if(coord[this->dimension-1]<=finish) goto start;
		// end of the loop

		//copy the data from vector result to array adjacentcelllist
		for(i=0; i<max; i++)
		{
			signed char * now=this->splitvoxellist+i*MaxDimension;
			int j;
			for(j=0; j<this->dimension; j++)
				now[j]=result[i].c[j];
			for(; j<MaxDimension; j++)
				now[j]=0;
		}

		//function end
		return;
}

VoxelList::VoxelList()
{
	this->dimension=0;
	this->level=0;
	this->halfsize=0;
	this->originhalfsize=0;
	this->nbrsplitvoxel=0;
	this->originindex=nullptr;
	this->subindex=nullptr;
	this->splitvoxellist=nullptr;
}
VoxelList::~VoxelList()
{
	if(this->originindex!=nullptr) delete this->originindex;
	if(this->subindex!=nullptr) delete this->subindex;
	if(this->splitvoxellist!=nullptr) delete [] this->splitvoxellist;
}
VoxelList::VoxelList(int Dimension, voxeltype VoxelRank, size_t ExpectedSize)
{
	this->dimension=Dimension;
	this->level=0;
	this->halfsize=::SystemSize/VoxelRank/2;
	this->originhalfsize=this->halfsize;
	this->getsplitvoxellist();
	this->originindex=new std::vector<voxeltype>;
	this->originindex->reserve(ExpectedSize*MaxDimension);
	this->subindex=new std::vector<unsigned short>;
}
void VoxelList::GetOriginalList(const SpherePacking & Config, double Radius)
{	
	double increment=2.0*this->halfsize;
	double finish=::SystemSize;
	double initial=this->halfsize;
	double coord[MaxDimension];
	short i;
	//initialize
	for(i=0;i<this->dimension;i++)
		coord[i]=initial;
	//main loop body
start:
	if(Config.CheckVoxelOverlap( Coord2Vect(coord), this->halfsize, Radius)==false)
	{
		for(i=0;i<this->dimension;i++)
			this->originindex->push_back((voxeltype)(coord[i]/increment));//convert voxel center coordinate back to original voxel index
		for(;i<MaxDimension;i++)
			this->originindex->push_back(0);
	}
	//loop end
	coord[0]+=increment;
	for(i=0;i<this->dimension-1;i++)
		if(coord[i]>finish)
		{
			coord[i]=initial;
			coord[i+1]+=increment;
		}
		if(coord[this->dimension-1]<finish) goto start;

		//function end
		return;
}
unsigned long VoxelList::NumVoxel(void) const
{
	return this->originindex->size()/MaxDimension;
}
void VoxelList::GetCenter(double * result, size_t nbr) const//get the center coordinate of the nbr-th voxel
{
	std::vector<voxeltype>::iterator itero=this->originindex->begin()+MaxDimension*nbr;
	std::vector<unsigned short>::iterator iters=this->subindex->begin()+this->level*nbr;
	double nowsize=this->originhalfsize;

	for(int j=0; j<this->dimension; j++)
	{
		result[j]=(2*itero[j]+1)*nowsize;//result is center of original voxel
		assert(result[j]<=1+1e-4);
		assert(result[j]>=-1e-4);
	}
	for(int i=0; i<this->level; i++)
	{
		signed char * iterc=this->splitvoxellist+iters[i]*MaxDimension;
		nowsize/=2;
		for(int j=0; j<this->dimension; j++)
		{
			result[j]+=nowsize*iterc[j];
			assert(result[j]<=1+1e-4);
			assert(result[j]>=-1e-4);
		}
	}//result is the center of subvoxel

}
void VoxelList::GetRandomCoordinate(double * result, RandomGenerator & gen) const//generate a random coordinate which is inside a random voxel, write it into result
{
	size_t nbr;
	do
	{
		nbr=(size_t)(gen.RandomDouble()*this->NumVoxel());
	}
	while(nbr>=this->NumVoxel());
	this->GetCenter(result, nbr);
	for(int j=0; j<this->dimension; j++)
	{
		double rand2=2*gen.RandomDouble()-1;
		result[j]+=rand2*this->halfsize;
		assert(result[j]<=1+1e-4);
		assert(result[j]>=-1e-4);
	}
}
bool VoxelList::CheckOverlapDeep(const SpherePacking & Config, const double & Radius, const double * Center, double QuaterSize, unsigned int checklevel)
{
	if(Config.CheckVoxelOverlap( Coord2Vect(Center), 2*QuaterSize, Radius)==true) return true;
	if(checklevel==0)
		return false;
	else
	{
		double newCenter[::MaxDimension];
		for(unsigned short i=0; i<this->nbrsplitvoxel; i++)
		{
			for(int j=0; j< ::MaxDimension; j++)
				newCenter[j]=Center[j]+QuaterSize*this->splitvoxellist[i*::MaxDimension+j];
			if(Config.CheckOverlap( Coord2Vect(newCenter), Radius)==false)
				return false;
			if(this->CheckOverlapDeep(Config, Radius, newCenter, QuaterSize/2, checklevel-1)==false)
				return false;
		}
		return true;
	}
}
void VoxelList::SplitVoxelSerial(const SpherePacking & Config, double Radius, int checklevel)
{
	std::vector<voxeltype> * neworigin = new std::vector<voxeltype>;
	std::vector<unsigned short> * newsub = new std::vector<unsigned short>;
	neworigin->reserve(this->NumVoxel()*MaxDimension);
	newsub->reserve(this->NumVoxel()*(this->level+1));
	unsigned long nbr=this->NumVoxel();

	for(signed long i=0; i<nbr; i++)
	{
		//locate this voxel
		std::vector<voxeltype>::iterator itero=this->originindex->begin()+MaxDimension*i;
		std::vector<unsigned short>::iterator iters=this->subindex->begin()+this->level*i;
		double center[MaxDimension], ncenter[MaxDimension];
		this->GetCenter(center, i);
		if(Config.CheckVoxelOverlap( Coord2Vect(center), this->halfsize, Radius)==true)
			continue;
		for(unsigned short j=0; j<this->nbrsplitvoxel; j++)
		{
			signed char * now=this->splitvoxellist+j*MaxDimension;
			//calculate new center
			for(int k=0; k<this->dimension; k++)
				ncenter[k]=center[k]+now[k]*this->halfsize*0.5;
			if(this->CheckOverlapDeep(Config, Radius, ncenter, this->halfsize/4, checklevel)==false)
			{
				//add new voxel

				int k;
				for(k=0; k<this->dimension; k++)
					neworigin->push_back(itero[k]);
				for(; k<MaxDimension; k++)
					neworigin->push_back(0);
				for(k=0; k<this->level; k++)
					newsub->push_back(iters[k]);

				newsub->push_back(j);

			}
		}
	}

	delete this->originindex;
	delete this->subindex;
	this->originindex=neworigin;
	this->subindex=newsub;
	this->halfsize/=2;
	this->level++;

}
void VoxelList::DistillVoxelList(const SpherePacking & Config, double Radius, int checklevel)
{
	std::vector<voxeltype> * neworigin = new std::vector<voxeltype>;
	std::vector<unsigned short> * newsub = new std::vector<unsigned short>;
	neworigin->reserve(this->NumVoxel()*MaxDimension);
	newsub->reserve(this->NumVoxel()*(this->level+1));

	omp_lock_t lock;
	omp_init_lock(&lock);
	signed long max=this->NumVoxel();
#pragma omp parallel default(shared)
	{
		SpherePacking myList(Config);
#pragma omp for
		for(signed long i=0; i<max; i++)
		{
			double center[MaxDimension];
			this->GetCenter(center, i);
			if(this->CheckOverlapDeep(myList, Radius, center, this->halfsize/2, checklevel)==true)
				continue;
			::omp_set_lock(&lock);
			for(int j=0; j< ::MaxDimension; j++)
				neworigin->push_back((*this->originindex)[i*::MaxDimension+j]);
			for(int j=0; j<this->level; j++)
				newsub->push_back((*this->subindex)[i*this->level+j]);
			::omp_unset_lock(&lock);
		}
	}
	::omp_destroy_lock(&lock);
	delete this->originindex;
	delete this->subindex;
	this->originindex=neworigin;
	this->subindex=newsub;
}
void VoxelList::SplitVoxel(const SpherePacking & Config, double Radius, int checklevel)
{
	std::vector<voxeltype> * neworigin = new std::vector<voxeltype>;
	std::vector<unsigned short> * newsub = new std::vector<unsigned short>;
	neworigin->reserve(this->NumVoxel()*MaxDimension);
	newsub->reserve(this->NumVoxel()*(this->level+1));
	unsigned long nbr=this->NumVoxel();

	omp_lock_t lock;
	//omp_init_lock(&lock);
#pragma omp parallel default(shared)
	{
		SpherePacking myList(Config);
		std::vector<voxeltype> myneworigin;
		std::vector<unsigned short> mynewsub;
#pragma omp for
		for(signed long i=0; i<nbr; i++)
		{
			//locate this voxel
			std::vector<voxeltype>::iterator itero=this->originindex->begin()+MaxDimension*i;
			std::vector<unsigned short>::iterator iters=this->subindex->begin()+this->level*i;
			double center[MaxDimension], ncenter[MaxDimension];
			this->GetCenter(center, i);
			if(myList.CheckVoxelOverlap( Coord2Vect(center), this->halfsize, Radius)==true)
				continue;
			for(unsigned short j=0; j<this->nbrsplitvoxel; j++)
			{
				signed char * now=this->splitvoxellist+j*MaxDimension;
				//calculate new center
				for(int k=0; k<this->dimension; k++)
					ncenter[k]=center[k]+now[k]*this->halfsize*0.5;
				if(this->CheckOverlapDeep(myList, Radius, ncenter, this->halfsize/4, checklevel)==false)
				{
					//add new voxel
					//omp_set_lock(&lock);
					int k;
					for(k=0; k<this->dimension; k++)
						myneworigin.push_back(itero[k]);
					for(; k<MaxDimension; k++)
						myneworigin.push_back(0);
					for(k=0; k<this->level; k++)
						mynewsub.push_back(iters[k]);

					mynewsub.push_back(j);
					//omp_unset_lock(&lock);
				}
			}
		}
#pragma omp critical
		{
			if(myneworigin.size()>0)
			{
				auto size=neworigin->size();
				neworigin->resize(size+myneworigin.size());
				std::memcpy(&(*neworigin)[size], &myneworigin[0], sizeof(voxeltype)*myneworigin.size());
			}
			if(mynewsub.size()>0)
			{
				auto size2=newsub->size();
				newsub->resize(size2+mynewsub.size());
				std::memcpy(&(*newsub)[size2], &mynewsub[0], sizeof(unsigned short)*mynewsub.size());
			}
		}
	}
	//omp_destroy_lock(&lock);

	delete this->originindex;
	delete this->subindex;
	this->originindex=neworigin;
	this->subindex=newsub;
	this->halfsize/=2;
	this->level++;

}
double VoxelList::GetVoxelSize() const
{
	return this->halfsize*2;
}



//convert double * to GeometryVector
GeometryVector RSA::Coord2Vect(const double * coord)
{
	GeometryVector cr(this->dimension);
	std::memcpy(cr.x, coord, this->dimension*sizeof(double));
	return cr;
}

//initialize square box with side length 1, sphere radii are set so that around NumSphere pSpheres exist in saturation limit
RSA::RSA(int Dimension, unsigned long NumSphere)
{
	this->dimension=Dimension;

	//calculate diameter
	double density=(0.1059*this->dimension*this->dimension+0.3571*this->dimension+1.039)/std::pow((double)(2),this->dimension);
	this->spherevolume=std::pow(::SystemSize,this->dimension)*density/NumSphere;
	//this->diameter=2*std::pow(this->spherevolume/ HyperSphere_Volume(this->dimension, 1.0),(double)(1)/this->dimension);
	this->radius=std::pow(this->spherevolume/ HyperSphere_Volume(this->dimension, 1.0),(double)(1)/this->dimension);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<Dimension; i++)
	{
		GeometryVector v(Dimension);
		v.x[i]=1.0;
		basis.push_back(v);
	}
	pSpheres=new SpherePacking(this->dimension, &basis[0], 2*this->radius);
	pVoxels=nullptr;
}

//initial configuration pa, upcoming pSpheres have one radii
//pa should support pSpheres of different sizes, but this is not tested.
RSA::RSA(const SpherePacking & pa, double radius)
{
	this->dimension=pa.GetDimension();

	this->radius=radius;
	this->spherevolume=HyperSphere_Volume(this->dimension, this->radius);

	pSpheres=new SpherePacking(pa);
	pVoxels=nullptr;
}
RSA::~RSA()
{
	if(this->pSpheres!=nullptr) delete this->pSpheres;
	if(this->pVoxels!=nullptr) delete this->pVoxels;
}

void RSA::RSA_I(size_t NumInsertedLimit, size_t TrialTime)//try to insert sphere, if less than NumInsertedLimit pSpheres inserted in TrialTime trials, stop
{
	size_t nbrinserted=10000;
	double temp[MaxDimension];
	while(nbrinserted>NumInsertedLimit)
	{
		nbrinserted=0;
		for(size_t i=0; i<TrialTime; i++)
		{
			for(int j=0; j<this->dimension; j++)
				temp[j]=this->gen.RandomDouble();
			if(this->pSpheres->CheckOverlap( Coord2Vect(temp), this->radius)==false)
			{
				this->pSpheres->Insert( this->radius, Coord2Vect(temp));
				nbrinserted++;
			}
		}
	}
}
void RSA::RSA_I(size_t NumCut)//try to insert sphere until there are NumCut spheres
{
	double temp[MaxDimension];
	while (this->pSpheres->NumParticle() < NumCut)
	{
		for (int j = 0; j < this->dimension; j++)
			temp[j] = this->gen.RandomDouble();
		if (this->pSpheres->CheckOverlap(Coord2Vect(temp), this->radius) == false)
		{
			this->pSpheres->Insert(this->radius, Coord2Vect(temp));
		}
	}
}
void RSA::RSA_II_Serial(size_t NumInsertedLimit, size_t TrialTime)
{
	if(this->pVoxels->NumVoxel()==0) return;
	size_t nbrinserted=10000;
	double temp[MaxDimension];
	while(nbrinserted>NumInsertedLimit)
	{
		nbrinserted=0;
		for(size_t i=0; i<TrialTime; i++)
		{
			for(int j=0; j<this->dimension; j++)
				this->pVoxels->GetRandomCoordinate(temp, this->gen);
			if(this->pSpheres->CheckOverlap( Coord2Vect(temp), this->radius)==false)
			{
				this->pSpheres->Insert( this->radius, Coord2Vect(temp));
				nbrinserted++;
			}
		}
	}
}	
struct RSA_II_TempStruct
{
	double x[MaxDimension];
};
class RSA_II_TempClass
{
public:
	RandomGenerator * g;
	RSA_II_TempClass(RandomGenerator & gen)
	{
		this->g= & gen;
	}
	ptrdiff_t operator() (ptrdiff_t max)
	{
		return static_cast<ptrdiff_t>(this->g->RandomDouble()*max);
	}
};
void RSA::RSA_II_Parallel( size_t NumInsertedLimit, size_t TrialTime)
{
	if(this->pVoxels->NumVoxel()==0) return;
	size_t nbrinserted=10000;
	std::vector<RSA_II_TempStruct> results;

	while(nbrinserted>NumInsertedLimit)
	{
		nbrinserted=0;
#pragma omp parallel default(shared)
		{
			SpherePacking myList(*this->pSpheres);
			RandomGenerator mygen((int)(this->gen.RandomDouble()*10000)+omp_get_thread_num());
			std::vector<RSA_II_TempStruct> myresults;
#pragma omp for
			for(signed long i=0; i<TrialTime; i++)
			{
				double temp[MaxDimension];
				for(int j=0; j<this->dimension; j++)
					this->pVoxels->GetRandomCoordinate(temp, mygen);
				if(myList.CheckOverlap(Coord2Vect(temp), this->radius)==false)
				{
					myresults.push_back(RSA_II_TempStruct());
					for(int i=0; i< ::MaxDimension; i++)
						myresults.back().x[i]=temp[i];

				}
			}

#pragma omp critical
			{
				while(myresults.size()!=0)
				{
					results.push_back(myresults.back());
					myresults.pop_back();
				}
			}
		}
		RSA_II_TempClass rtemp(this->gen);
		std::random_shuffle(results.begin(), results.end(), rtemp);
		for(auto iter=results.begin(); iter!=results.end(); iter++)
		{
			double * temp = iter->x;
			if(this->pSpheres->CheckOverlap(Coord2Vect(temp), this->radius)==false)
			{
				this->pSpheres->Insert(this->radius, Coord2Vect(temp));
				nbrinserted++;
			}
		}

		//std::cout<<"size of temp rsults:"<<results.size()<<'\n';

		results.clear();
	}
}	
void RSA::RSA_II(size_t NumInsertedLimit, size_t TrialTime, bool parallel)
{
	if(parallel)
		this->RSA_II_Parallel(NumInsertedLimit, TrialTime);
	else
		this->RSA_II_Serial(NumInsertedLimit, TrialTime);
}
double RSA::Density(void) const
{
	return this->pSpheres->NumParticle() * ::HyperSphere_Volume(this->dimension, this->radius);
}
void RSA::GetVoxelList(double VoxelDiameterRatio)//get original voxel list
{
	unsigned long rank=std::floor( std::sqrt(this->pSpheres->GetBasisVector(0).Modulus2())/this->radius/2.0/VoxelDiameterRatio);
	this->pVoxels=new VoxelList(this->dimension, rank);
	this->pVoxels->GetOriginalList(* this->pSpheres, this->radius);
}
void RSA::SplitVoxel(bool parallel, int checklevel)
{
	if(parallel)
		this->pVoxels->SplitVoxel(*this->pSpheres, this->radius, checklevel);
	else
		this->pVoxels->SplitVoxelSerial(*this->pSpheres, this->radius, checklevel);
}
void RSA::OutputSpheres(std::ostream & out)
{
	out<<this->radius;
	::Output(out, Configuration(*this->pSpheres, "A"));
}


SpherePacking GenerateRSAPacking(int dimension, unsigned long nbrsphere, short nbrinsertedlimit, unsigned long trial1, unsigned long trial2, double voxelratio, int seed, std::ostream & output, double * volumeratio, bool Verbose, size_t NumCut)
{
	unsigned long OriginalTrial2=trial2;
	bool tryagain=false;
	bool parallel=false;
	bool printconfig=true;
	if(dimension<0)//use dimension < 0 to indicate supress printing the configuration
	{
		dimension*=-1;
		printconfig=false;
	}
	if(nbrinsertedlimit<0)//use nbrinsertedlimit<0 to indicate parallel
	{
		nbrinsertedlimit*=-1;
		parallel=true;
	}
	if(voxelratio<0)//use voxelratio<0 to indicate do a voxel number test at saturation limit
	{
		voxelratio*=-1;
		tryagain=true;
	}


	RSA con1(dimension, nbrsphere);
	con1.gen.seed(seed);
	
	con1.RSA_I(nbrinsertedlimit, trial1);
	output<<"at time"<<std::time(NULL)- ::ProgramStart<<", "<<con1.pSpheres->NumParticle()<<" Spheres inserted\n";

	if(con1.pSpheres->NumParticle() >= NumCut)
	{
		while(con1.pSpheres->NumParticle() > NumCut)
			con1.pSpheres->DeleteParticle(con1.pSpheres->NumParticle()-1);
		output<<"Return unsaturated packing with "<<con1.pSpheres->NumParticle()<<" particles.\n";
		return SpherePacking(*con1.pSpheres);
	}

	con1.GetVoxelList(voxelratio);
	if(Verbose)
	{
		output<<"voxel size:"<<con1.pVoxels->GetVoxelSize()<<" \tvoxel centers are:\n";
		for(size_t i=0; i<con1.pVoxels->NumVoxel(); i++)
		{
			double temp[MaxDimension];
			con1.pVoxels->GetCenter(temp, i);
			for(size_t j=0; j<dimension; j++)
				output<<temp[j]<<", \t";
			output<<",\n";
		}
	}

	int SplitCount=0;
	while(con1.pVoxels->NumVoxel()!=0)
	{
		output<<"at time"<<std::time(NULL)- ::ProgramStart<<", "<<"current voxel:"<<con1.pVoxels->NumVoxel()<<'\n';
		con1.RSA_II(nbrinsertedlimit, trial2, parallel);
		output<<"at time"<<std::time(NULL)- ::ProgramStart<<", "<<con1.pSpheres->NumParticle()<<" Spheres inserted\n";

		if(con1.pSpheres->NumParticle() >= NumCut)
		{
			while(con1.pSpheres->NumParticle() > NumCut)
				con1.pSpheres->DeleteParticle(con1.pSpheres->NumParticle()-1);
			output<<"Return unsaturated packing with "<<con1.pSpheres->NumParticle()<<" particles.\n";
			return SpherePacking(*con1.pSpheres);
		}

		unsigned long temp=con1.pVoxels->NumVoxel();
		con1.SplitVoxel(parallel, 1);
		double ratio=static_cast<double>(con1.pVoxels->NumVoxel())/temp;
		if(ratio<1)
			trial2*=ratio;
		if(ratio>1.5)
			trial2=OriginalTrial2;
		SplitCount++;

		if(Verbose)
		{
			output<<"voxel size:"<<con1.pVoxels->GetVoxelSize()<<" \tvoxel centers are:\n";
			for(size_t i=0; i<con1.pVoxels->NumVoxel(); i++)
			{
				double temp[MaxDimension];
				con1.pVoxels->GetCenter(temp, i);
				for(size_t j=0; j<dimension; j++)
					output<<temp[j]<<", \t";
				output<<",\n";
			}
		}
	}
	if(tryagain)
	{
		delete con1.pVoxels;
		con1.GetVoxelList(voxelratio);
		while(con1.pVoxels->NumVoxel()!=0)
		{
			output<<"at time"<<std::time(NULL)- ::ProgramStart<<", "<<"current voxel:"<<con1.pVoxels->NumVoxel()<<'\n';
			con1.SplitVoxel(parallel, 1);
		}
	}


	output<<"at time"<<std::time(NULL)- ::ProgramStart<<", "<<con1.pSpheres->NumParticle()<<" Spheres inserted\n";
	output<<"r="<<con1.radius<<"\nVolume of a Sphere:"<<con1.spherevolume<<"\nDensity:"<<con1.Density()<<'\n';
	if(volumeratio!=nullptr)
		*volumeratio=con1.spherevolume;

	if(printconfig)
	{
		std::fstream outfile("SphereCoordinate.txt", std::ios::out);
		outfile.precision(17);
		con1.OutputSpheres(outfile);
	}

	return SpherePacking(*con1.pSpheres);
}
