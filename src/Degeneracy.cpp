#include "Degeneracy.h"
#include "etc.h"
#include "Solvers.h"
#include "ThreeBody.h"
#include <omp.h>


////////////////////////////////////////////////////////////////////////////
//codes related to searching degenerate configurations
#include <nlopt.h>
struct Aux2
{
	size_t NumParticle1;
	size_t NumParticle2;
	DimensionType Dim;
	double *pMinDisntace;
	size_t NumEvaluation;
};
Configuration GetConfig1(Aux2 * pAux, const double * x)
{
	DimensionType & dim=pAux->Dim;
	assert(dim< ::MaxDimension);
	GeometryVector base[::MaxDimension];
	for(DimensionType i=0; i<dim; i++)
	{
		base[i]=GeometryVector(dim);
		for(DimensionType j=0; j<dim; j++)
			base[i].x[j]=x[i*dim+j];
	}
	Configuration list(dim, base, 1.0);
	for(size_t i=0; i<pAux->NumParticle1; i++)
	{
		GeometryVector ParticleRelative(dim);
		for(DimensionType j=0; j<dim; j++)
			ParticleRelative.x[j]=x[dim*dim+dim*i+j];
		list.Insert("H", ParticleRelative);
	}
	return list;
}
Configuration GetConfig2(Aux2 * pAux, const double * x)
{
	DimensionType & dim=pAux->Dim;
	x=x+dim*dim+dim*pAux->NumParticle1;
	assert(dim< ::MaxDimension);
	GeometryVector base[::MaxDimension];
	for(DimensionType i=0; i<dim; i++)
	{
		base[i]=GeometryVector(dim);
		for(DimensionType j=0; j<dim; j++)
			base[i].x[j]=x[i*dim+j];
	}
	Configuration list(dim, base, 1.0);
	for(size_t i=0; i<pAux->NumParticle2; i++)
	{
		GeometryVector ParticleRelative(dim);
		for(DimensionType j=0; j<dim; j++)
			ParticleRelative.x[j]=x[dim*dim+dim*i+j];
		list.Insert("H", ParticleRelative);
	}
	return list;
}

double Objective2(unsigned n, const double *x, double *grad, void *data)
{
	Aux2 * pAux=reinterpret_cast<Aux2 *>(data);
	DimensionType & dim=pAux->Dim;

	Configuration list1 = GetConfig1(pAux, x);
	for(DimensionType i=0; i<dim; i++)
	{
		GeometryVector t=list1.GetBasisVector(i);
		if(t.Modulus2()<(*pAux->pMinDisntace)*(*pAux->pMinDisntace))
			return 10.0;
	}

	Configuration list2 = GetConfig2(pAux, x);
	for(DimensionType i=0; i<dim; i++)
	{
		GeometryVector t=list2.GetBasisVector(i);
		if(t.Modulus2()<(*pAux->pMinDisntace)*(*pAux->pMinDisntace))
			return 10.0;
	}

	assert(grad == nullptr);
	pAux->NumEvaluation++;

	FromStructure s1(list1);
	FromStructure s2(list2);
	double result=SolverDistance(&s1, &s2);
	//debug temp
	//std::cout<<"n="<<pAux->NumEvaluation<<"f="<<result<<'\n';
	return result;
}

void OptimizeDegeneracy2(Configuration & List1, Configuration & List2, double MinDistance)
{
	DimensionType dim=List1.GetDimension();
	size_t Num1=dim*dim+dim*List1.NumParticle();
	size_t NumParam=Num1+dim*dim+dim*List2.NumParticle();
	double * pParams = new double[NumParam];
	double * LBounds = new double[NumParam];
	double * UBounds = new double[NumParam];
	double * StepSizes = new double[NumParam];
	//double * Tolerences = new double[NumParam];

	//set basis vector parameters
	for(DimensionType i=0; i<dim; i++)
	{
		GeometryVector vBase = List1.GetBasisVector(i);
		double Length = std::sqrt(vBase.Modulus2());
		for(DimensionType j=0; j<dim; j++)
		{
			pParams[i*dim+j]=vBase.x[j];
			LBounds[i*dim+j]=vBase.x[j]-Length;
			UBounds[i*dim+j]=vBase.x[j]+Length;
			StepSizes[i*dim+j]=0.001*Length;
			//Tolerences[i*dim+j]=1e-10*Length;
		}
	}
	for(size_t i=0; i<List1.NumParticle(); i++)
	{
		for(DimensionType j=0; j<dim; j++)
		{
			//Configuration::particle * pA=List1.GetParticle(i);
			//pParams[dim*dim+i*dim+j] = pA->RelativeCoordinates.x[j];
			pParams[dim*dim+i*dim+j] = List1.GetRelativeCoordinates(i).x[j];
			LBounds[dim*dim+i*dim+j] = 0.0;
			UBounds[dim*dim+i*dim+j] = 1.0;
			StepSizes[dim*dim+i*dim+j] = 0.001;
			//Tolerences[dim*dim+i*dim+j]=1e-10;
		}
	}
	for(DimensionType i=0; i<dim; i++)
	{
		GeometryVector vBase = List2.GetBasisVector(i);
		double Length = std::sqrt(vBase.Modulus2());
		for(DimensionType j=0; j<dim; j++)
		{
			pParams[Num1+i*dim+j]=vBase.x[j];
			LBounds[Num1+i*dim+j]=vBase.x[j]-Length;
			UBounds[Num1+i*dim+j]=vBase.x[j]+Length;
			StepSizes[Num1+i*dim+j]=0.001*Length;
			//Tolerences[Num1+i*dim+j]=1e-10*Length;
		}
	}
	for(size_t i=0; i<List2.NumParticle(); i++)
	{
		for(DimensionType j=0; j<dim; j++)
		{
			pParams[Num1+dim*dim+i*dim+j] = List2.GetRelativeCoordinates(i).x[j];
			LBounds[Num1+dim*dim+i*dim+j] = 0.0;
			UBounds[Num1+dim*dim+i*dim+j] = 1.0;
			StepSizes[Num1+dim*dim+i*dim+j] = 0.001;
			//Tolerences[Num1+dim*dim+i*dim+j]=1e-10;
		}
	}

	//do optimization
	Aux2 aux;
	aux.Dim=dim;
	aux.NumParticle1=List1.NumParticle();
	aux.NumParticle2=List2.NumParticle();
	aux.pMinDisntace= & MinDistance;
	aux.NumEvaluation=0;
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LN_NEWUOA, NumParam);

	nlopt_set_lower_bounds(opt, LBounds);
	nlopt_set_upper_bounds(opt, UBounds);
	nlopt_set_initial_step(opt, StepSizes);
	nlopt_set_min_objective(opt, Objective2, reinterpret_cast<void *>(&aux));
	//nlopt_set_xtol_abs(opt, Tolerences);
	//nlopt_set_ftol_rel(opt, 1e-10);
	nlopt_set_xtol_rel(opt, 1e-15);
	nlopt_set_maxtime(opt, 300.0);

	double minf; /* the minimum objective value, upon return */

	//std::cout<<"NLopt fvalue before:"<<Objective(NumParam, pParams, nullptr, reinterpret_cast<void *>(&aux))<<'\n';

	auto status=nlopt_optimize(opt, pParams, &minf);

	//debug temp
	//std::cout<<"NLopt return:"<<status<<'\n';
	//std::cout<<"NLopt fvalue after:"<<minf<<'\n';

	List1 = GetConfig1(&aux, pParams);
	List2 = GetConfig2(&aux, pParams);
	nlopt_destroy(opt);
	delete [] pParams;
	delete [] LBounds;
	delete [] UBounds;
	delete [] StepSizes;
	//delete [] Tolerences;
}
void SearchDegenerate2(RandomGenerator & gen, size_t SearchTime, std::vector<std::pair<Configuration, Configuration> > & results, size_t NumParticle1, size_t NumParticle2, DimensionType dim, double * pd2, double * pd3)
{
	for(signed long nn=0; nn<SearchTime; nn++)
	{
		//initialize a configuration
		unsigned short cellrank = std::floor( std::pow(NumParticle1/ParticlePerCell, static_cast<double>(1)/dim) );//about 1 particles per cell
		double initsize=1.0;
		std::vector<GeometryVector> base;
		for(DimensionType i=0; i<dim; i++)
		{
			base.push_back(GeometryVector(dim));
			base.back().x[i]=initsize;
			for(DimensionType j=i+1; j<dim; j++)
				base.back().x[j]=initsize*(gen.RandomDouble()-0.5);
		}
		Configuration conf1(dim, &base[0], initsize/cellrank);
		for(size_t i=0; i<NumParticle1; i++)
		{
			GeometryVector temp;
			for(DimensionType j=0; j<dim; j++)
				temp.x[j]=gen.RandomDouble();
			conf1.Insert("C", temp);
		}

		base.clear();
		for(DimensionType i=0; i<dim; i++)
		{
			base.push_back(GeometryVector(dim));
			base.back().x[i]=initsize;
			for(DimensionType j=i+1; j<dim; j++)
				base.back().x[j]=initsize*(gen.RandomDouble()-0.5);
		}
		Configuration conf2(dim, &base[0], initsize/cellrank);
		for(size_t i=0; i<NumParticle2; i++)
		{
			GeometryVector temp;
			for(DimensionType j=0; j<dim; j++)
				temp.x[j]=gen.RandomDouble();
			conf2.Insert("C", temp);
		}

		OptimizeDegeneracy2(conf1, conf2, 0.0);
		FromStructure Optimized1(conf1);
		FromStructure Optimized2(conf2);
		double d2=SolverDistance(&Optimized1, &Optimized2);
		if (pd2 != nullptr)
			(*pd2) = d2;

		//std::cout<<"dist="<<TwoBodyDistance<<" \t";
		if(d2<1e-7)
		{
			//two body degeneracy found
			double d3=ThreeBodySolverDistance(&Optimized1, &Optimized2);
			if (pd3 != nullptr)
				(*pd3) = d3;

			std::cout<<"Found a degenerate structure, 2-body distance:"<<d2<<", 3-body distance:"<<d3<<'\n';

			results.push_back(std::pair<Configuration, Configuration>(conf1, conf2));
		}
		//std::cout<<'\n';
	}
}











struct Aux
{
	size_t NumParticle;
	DimensionType Dim;
	double *pMinDisntace;
	size_t NumEvaluation;
	LatticeSumSolver * pReference;
};
Configuration GetConfig(Aux * pAux, const double * x)
{
	DimensionType & dim=pAux->Dim;
	assert(dim< ::MaxDimension);
	GeometryVector base[::MaxDimension];
	for(DimensionType i=0; i<dim; i++)
	{
		base[i]=GeometryVector(dim);
		for(DimensionType j=0; j<dim; j++)
			base[i].x[j]=x[i*dim+j];
	}
	Configuration list(dim, base, 1.0);
	for(size_t i=0; i<pAux->NumParticle; i++)
	{
		GeometryVector ParticleRelative(dim);
		for(DimensionType j=0; j<dim; j++)
			ParticleRelative.x[j]=x[dim*dim+dim*i+j];
		list.Insert("H", ParticleRelative);
	}
	return list;
}

double Objective(unsigned n, const double *x, double *grad, void *data)
{
	Aux * pAux=reinterpret_cast<Aux *>(data);
	DimensionType & dim=pAux->Dim;
	assert(n==dim*dim+dim*pAux->NumParticle);

	Configuration list = GetConfig(pAux, x);
	double Volume = list.PeriodicVolume();
	if(Volume< std::pow(::LengthPrecision, static_cast<double>(dim)))
		return ::MaxEnergy;

	for(DimensionType i=0; i<dim; i++)
	{
		GeometryVector t=list.GetBasisVector(i);
		if(t.Modulus2()<(*pAux->pMinDisntace)*(*pAux->pMinDisntace))
			return 10.0;
	}

	assert(grad == nullptr);
	pAux->NumEvaluation++;

	FromStructure s(list);
	double result=SolverDistance(pAux->pReference, &s);
	//debug temp
	//std::cout<<"n="<<pAux->NumEvaluation<<"f="<<result<<'\n';
	return result;
}

void OptimizeDegeneracy(Configuration & List, LatticeSumSolver * pReference, double MinDistance)
{
	DimensionType dim=List.GetDimension();
	size_t NumParam=dim*dim+dim*List.NumParticle();
	double * pParams = new double[NumParam];
	double * LBounds = new double[NumParam];
	double * UBounds = new double[NumParam];
	double * StepSizes = new double[NumParam];
	//double * Tolerences = new double[NumParam];

	//set basis vector parameters
	for(DimensionType i=0; i<dim; i++)
	{
		GeometryVector vBase = List.GetBasisVector(i);
		double Length = std::sqrt(vBase.Modulus2());
		for(DimensionType j=0; j<dim; j++)
		{
			pParams[i*dim+j]=vBase.x[j];
			LBounds[i*dim+j]=vBase.x[j]-0.5*Length;
			UBounds[i*dim+j]=vBase.x[j]+0.5*Length;
			StepSizes[i*dim+j]=0.001*Length;
			//Tolerences[i*dim+j]=1e-10*Length;
		}
	}
	for(size_t i=0; i<List.NumParticle(); i++)
	{
		for(DimensionType j=0; j<dim; j++)
		{
			pParams[dim*dim+i*dim+j] = List.GetRelativeCoordinates(i).x[j];
			LBounds[dim*dim+i*dim+j] = 0.0;
			UBounds[dim*dim+i*dim+j] = 1.0;
			StepSizes[dim*dim+i*dim+j] = 0.001;
			//Tolerences[dim*dim+i*dim+j]=1e-10;
		}
	}

	//do optimization
	Aux aux;
	aux.Dim=dim;
	aux.NumParticle=List.NumParticle();
	aux.pMinDisntace= & MinDistance;
	aux.NumEvaluation=0;
	aux.pReference=pReference;
	nlopt_opt opt;

	opt = nlopt_create(NLOPT_LN_NEWUOA, NumParam);

	nlopt_set_lower_bounds(opt, LBounds);
	nlopt_set_upper_bounds(opt, UBounds);
	nlopt_set_initial_step(opt, StepSizes);
	nlopt_set_min_objective(opt, Objective, reinterpret_cast<void *>(&aux));
	//nlopt_set_xtol_abs(opt, Tolerences);
	//nlopt_set_ftol_rel(opt, 1e-10);
	nlopt_set_xtol_rel(opt, 1e-15);
	nlopt_set_maxtime(opt, 300.0);

	double minf; /* the minimum objective value, upon return */

	//std::cout<<"NLopt fvalue before:"<<Objective(NumParam, pParams, nullptr, reinterpret_cast<void *>(&aux))<<'\n';

	auto status=nlopt_optimize(opt, pParams, &minf);

	//debug temp
	//std::cout<<"NLopt return:"<<status<<'\n';
	//std::cout<<"NLopt fvalue after:"<<minf<<'\n';

	List = GetConfig(&aux, pParams);
	nlopt_destroy(opt);
	delete [] pParams;
	delete [] LBounds;
	delete [] UBounds;
	delete [] StepSizes;
	//delete [] Tolerences;
}
void SearchDegenerate(LatticeSumSolver * pReference, RandomGenerator & GlobalGen, size_t SearchTime, std::vector<Configuration> & results, size_t NumParticle)
{
	static size_t Debug=0;
	//pre-allocate theta series of pReference
	SolverDistance(pReference, pReference);

	Configuration ref=pReference->GetStructure();
	DimensionType dim=ref.GetDimension();
	//size_t NumParticle=ref.NumParticle();

	//generate some terms
	double tempDist=1;
	while(pReference->Terms.size()<2)
	{
		pReference->UpdateTerms(tempDist);
		tempDist*=2;
	}

	double MinDistance=pReference->Terms[1].distance/2.0;
	if(MinDistance>0.5*std::sqrt(ref.GetBasisVector(0).Modulus2()))
		MinDistance=0.5*std::sqrt(ref.GetBasisVector(0).Modulus2());

	std::vector<int> seeds;
	for(signed long nn=0; nn<SearchTime; nn++)
		seeds.push_back(GlobalGen.RandomDouble()*10000);

#pragma omp parallel for schedule(dynamic) 
	for(signed long nn=0; nn<SearchTime; nn++)
	{
		RandomGenerator gen(seeds[nn]);
		//initialize a configuration
		std::cout<<'-';
		unsigned short cellrank = std::floor( std::pow(NumParticle/ParticlePerCell, static_cast<double>(1)/dim) );//about 1 particles per cell
		double initsize=std::pow(ref.PeriodicVolume(), static_cast<double>(1)/dim);
		std::vector<GeometryVector> base;
		for(DimensionType i=0; i<dim; i++)
		{
			base.push_back(GeometryVector(dim));
			base.back().x[i]=initsize;
			for(DimensionType j=i+1; j<dim; j++)
				base.back().x[j]=initsize*(gen.RandomDouble()-0.5);
		}
		Configuration conf(dim, &base[0], initsize/cellrank);
		size_t trialcount=0;
		for(size_t i=0; i<NumParticle;)
		{
			GeometryVector temp;
			for(DimensionType j=0; j<dim; j++)
				temp.x[j]=gen.RandomDouble();
			bool TooNearNeighbor=false;
			conf.IterateThroughNeighbors(temp, MinDistance, [&TooNearNeighbor](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
			{
				TooNearNeighbor=true;
			}, &TooNearNeighbor);
			if(TooNearNeighbor==false)
			{
				conf.Insert("C", temp);
				i++;
			}
			trialcount++;
			if(trialcount%100000==99999)
			{
				GeometryVector tbase[ ::MaxDimension];
				for(DimensionType j=0; j<dim; j++)
					tbase[j]=conf.GetBasisVector(j)*2;
				conf.ChangeBasisVector(tbase);
			}
		}

		OptimizeDegeneracy(conf, pReference, MinDistance);
		FromStructure Optimized(conf);
#pragma omp critical
		{
			if(SolverDistance(pReference, &Optimized)<1e-7)
			{
				//std::string name("Debug");
				//name+=std::to_string((long double)(Debug++));
				//::Plot(name, MultiplicateStructure(conf, 3));
				//::Output(name, conf);
				//two body degeneracy found
				if(ThreeBodySolverDistance(pReference, &Optimized)<1e-5)
				{
					std::cout<<"Found the target!\n";
				}
				else
				{
					bool HasFoundPreviously=false;
					for(auto iter=results.begin(); iter!=results.end(); iter++)
					{
						FromStructure ts(*iter);
						if(ThreeBodySolverDistance(&ts, &Optimized)<1e-5)
						{
							std::cout<<"Found some degenerate structure again!\n";
							HasFoundPreviously=true;
							break;
						}
					}
					if(HasFoundPreviously==false)
					{
						std::cout<<"Found a degenerate structure, 2-body distance:"<<SolverDistance(pReference, &Optimized)<<", 3-body distance:"<<ThreeBodySolverDistance(pReference, &Optimized)<<'\n';
						results.push_back(conf);
					}

				}
			}
		}
	}
}



Configuration FindStructure_inner(LatticeSumSolver * pReference, RandomGenerator & GlobalGen, size_t SearchTime, size_t NumParticle, double Precision)
{
	//pre-allocate theta series of pReference
	SolverDistance(pReference, pReference);

	Configuration ref=pReference->GetStructure();
	DimensionType dim=ref.GetDimension();
	//size_t NumParticle=ref.NumParticle();

	//generate some terms
	double tempDist=1;
	while(pReference->Terms.size()<2)
	{
		pReference->UpdateTerms(tempDist);
		tempDist*=2;
	}

	double MinDistance=pReference->Terms[1].distance/2.0;
	if(MinDistance>0.5*std::sqrt(ref.GetBasisVector(0).Modulus2()))
		MinDistance=0.5*std::sqrt(ref.GetBasisVector(0).Modulus2());

	std::vector<int> seeds;
	for(signed long nn=0; nn<SearchTime; nn++)
		seeds.push_back(GlobalGen.RandomDouble()*10000);

	volatile bool ResultFound=false;
	Configuration result=pReference->GetStructure();
	volatile double MinSolverDistance=100.0;
	size_t SameMinDistanceHit=0;

	progress_display pd(SearchTime);

#pragma omp parallel for schedule(dynamic) 
	for(signed long nn=0; nn<SearchTime; nn++)
	{
		if(ResultFound)
			continue;
		RandomGenerator gen(seeds[nn]);
		//initialize a configuration
		unsigned short cellrank = std::floor( std::pow(NumParticle/ParticlePerCell, static_cast<double>(1)/dim) );//about 1 particles per cell
		double initsize=std::pow(ref.PeriodicVolume(), static_cast<double>(1)/dim);
		std::vector<GeometryVector> base;
		for(DimensionType i=0; i<dim; i++)
		{
			base.push_back(GeometryVector(dim));
			base.back().x[i]=initsize;
			for(DimensionType j=i+1; j<dim; j++)
				base.back().x[j]=initsize*(gen.RandomDouble()-0.5);
		}
		Configuration conf(dim, &base[0], initsize/cellrank);
		size_t trialcount=0;
		for(size_t i=0; i<NumParticle;)
		{
			GeometryVector temp;
			for(DimensionType j=0; j<dim; j++)
				temp.x[j]=gen.RandomDouble();
			bool TooNearNeighbor=false;
			conf.IterateThroughNeighbors(temp, MinDistance, [&TooNearNeighbor](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
			{
				TooNearNeighbor=true;
			}, &TooNearNeighbor);
			if(TooNearNeighbor==false)
			{
				conf.Insert("C", temp);
				i++;
			}
			trialcount++;
			if(trialcount%100000==99999)
			{
				GeometryVector tbase[ ::MaxDimension];
				for(DimensionType j=0; j<dim; j++)
					tbase[j]=conf.GetBasisVector(j)*2;
				conf.ChangeBasisVector(tbase);
			}
		}

		OptimizeDegeneracy(conf, pReference, MinDistance);
		FromStructure Optimized(conf);
#pragma omp critical
		{
			pd++;
			double dist2=SolverDistance(pReference, &Optimized);
			if(MinSolverDistance>dist2)
			{
				MinSolverDistance=dist2;
				SameMinDistanceHit=1;
			}
			else if( (dist2-MinSolverDistance)/MinSolverDistance < 1e-6 )
				SameMinDistanceHit++;
			//debug temp
			logfile<<"d2="<<dist2<<'\n';
			if(dist2<Precision)
			{
				std::cout<<"find 2-body degeneracy\n";
				//double dist3=ThreeBodySolverDistance(pReference, &Optimized);
				//std::cout<<"d3="<<dist3<<'\n';
				double dist3=0;
				if(dist3<1e-3)
				{
				    std::cout<<"find 3-body degeneracy\n";
					result = conf;
					ResultFound=true;
				}
			}
		}
	}
	if(!ResultFound)
		std::cout<<'\n'<<NumParticle<<"-particle basis search did not find crystal structure. Minimum 2-body distance="<<MinSolverDistance<<" hit "<<SameMinDistanceHit<<" times\n";
	return result;
}

Configuration FindStructure(Configuration Reference, RandomGenerator & GlobalGen, size_t SearchTime, double Precision)
{
	double rho=Reference.NumParticle()/Reference.PeriodicVolume();
	Reference.Rescale( std::pow(rho, 1.0/Reference.GetDimension() ) );
	FromStructure sRef(Reference);
	size_t NumOriginal = Reference.NumParticle();
	for(int NumP=1; NumP<NumOriginal; NumP++)
	{
		Configuration c=FindStructure_inner(&sRef, GlobalGen, SearchTime, NumP, Precision);
		if(c.NumParticle()<NumOriginal)
			return c;
	}
	return Reference;
}
Configuration FindStructure_SpecificNumberBasis(Configuration Reference, RandomGenerator & GlobalGen, size_t NumBasis, size_t SearchTime, double Precision)
{
	double rho=Reference.NumParticle()/Reference.PeriodicVolume();
	Reference.Rescale( std::pow(rho, 1.0/Reference.GetDimension() ) );
	FromStructure sRef(Reference);
	size_t NumOriginal = Reference.NumParticle();
	Configuration c=FindStructure_inner(&sRef, GlobalGen, SearchTime, NumBasis, Precision);
	if(c.NumParticle()<NumOriginal)
		return c;
	return Reference;
}