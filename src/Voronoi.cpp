#include "Voronoi.h"
#include "DisplaySpheres.h"
#include "PairCorrelation.h"
#ifdef USE_MATHGL
#include <mgl2/mgl.h>
#endif

#ifdef USE_CGAL

#include <CGAL/Cartesian_d.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Kernel_d/Sphere_d.h>
#include <CGAL/Delaunay_d.h>

//#define USE_CGAL_PRECISE

typedef CGAL::Cartesian_d<CGAL::Gmpq> Kernele;
#ifdef USE_CGAL_PRECISE 
typedef CGAL::Homogeneous_d<CGAL::Gmpq> Kerneli;
#else
typedef CGAL::Homogeneous_d<double> Kerneli;
#endif
typedef CGAL::Sphere_d<Kerneli> Sphere;
typedef CGAL::Convex_hull_d<Kernele> ConvexHull;
typedef ConvexHull::Point_d cPoint;
typedef ConvexHull::Vertex_handle cVertex_handle;
typedef ConvexHull::Simplex_handle cSimplex_handle;
typedef CGAL::Delaunay_d<Kerneli> Delaunay;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Simplex_handle Simplex_handle;
typedef Delaunay::Point_d Point;

template <typename point> GeometryVector CGALPoint2GeometryVector(const point & p, DimensionType d)
{
	GeometryVector result(d);
	for(DimensionType i=0; i<d; i++)
		result.x[i]=CGAL::to_double(p.cartesian(i));
	return result;
}

// test if point "tested" is on the convex hull
bool IsOnConvexHull(const GeometryVector & tested, const std::vector<GeometryVector> & neighbors, DimensionType d)
{
	ConvexHull ch(d);
	cPoint p(d, tested.x, tested.x+d);
	cVertex_handle hTested = ch.insert(p);
	for(auto iter=neighbors.begin(); iter!=neighbors.end(); iter++)
		ch.insert(cPoint(d, iter->x, iter->x+d));
	for(ConvexHull::Hull_vertex_iterator iter=ch.hull_vertices_begin(); iter!=ch.hull_vertices_end(); iter++)
		if(iter->point() == hTested->point())
			return true;
	return false;
}

//calculate all vertexes of the voronoi cell of particle i, in vertices
//if pLinkage is not nullptr, it will be filled with linkage information of the format:
// (*pLinkage)[i*(Dimension+1)+j] is the index of a vertex linked with vertex i
// (*pLinkage)[i*(Dimension+1)+j] is -1 if vertex i have less than j linked vertexes.
void GetVoronoiCell(const PeriodicCellList<Empty> & c, size_t i, std::vector<GeometryVector> & vertices, std::vector<signed long> * pLinkage, std::set<size_t> * pNeighbors)
{
	DimensionType d=c.GetDimension();
	double radius=std::pow(c.PeriodicVolume()/c.NumParticle(), 1.0/d)*1.8;
	double FurthestVertexDist = ::MaxDistance;
	GeometryVector ci=c.GetCartesianCoordinates(i);
	GeometryVector ri=c.GetRelativeCoordinates(i);
	struct VoronoiVertexInfo
	{
		Simplex_handle hDelaunaySimplex;
		size_t VerticesVectorLocation;
	};
	std::vector<VoronoiVertexInfo> vInfo;
	vertices.clear();

	//debug temp
	//size_t NumIteration = 0;
	while(2*FurthestVertexDist>radius || vertices.empty())
	{
		//debug temp
		//NumIteration++;

		vInfo.clear();
		vertices.clear();
		radius*=1.25;
		if (pNeighbors != nullptr)
			pNeighbors->clear();
		//Get neighbors within this radius
		std::vector<GeometryVector> neighbors;
		c.IterateThroughNeighbors(ri, radius, [&i, &ci, &neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
		{
			if(Sourceparticle==i && shift.Modulus2()==0)
				return;
			neighbors.push_back(ci+shift);
		});
		//if no neighbors are found, then obviously radius is too small
		if(neighbors.empty())
			continue;
		//convex hull test, if the center particle is on convex hull then we are sure radius is not enough.
		if(IsOnConvexHull(ci, neighbors, d))
			continue;
		//Get Delaunay of neighbors
		Delaunay De(d);
		Point p(d, ci.x, ci.x+d);
		Vertex_handle hvi=De.insert(p);
		for(auto iter=neighbors.begin(); iter!=neighbors.end(); iter++)
		{
			Point p(d, iter->x, iter->x+d);
			De.insert(p);
		}
		//for each simplex in Delaunay, if one of its vertex is hvi, then calculate the circumsphere center.
		std::list<Simplex_handle> simplices = De.all_simplices();
		for(auto iter=simplices.begin(); iter!=simplices.end(); iter++)
		{
			bool has_hvi = false;
			bool IsValid = true;
			std::vector<Point> points;
			for(size_t i=0; i<=d; i++)
			{
				Vertex_handle hv=De.vertex_of_simplex(*iter, i);
				if(hv==hvi)
					has_hvi=true;
				if(hv!=nullptr)
					points.push_back(De.associated_point(hv));
				else
					IsValid=false;
			}
			if(IsValid==false)
				continue;
			if(has_hvi)
			{
				//calculate sphere center
				Sphere s(d, points.begin(), points.end());
				if(s.is_legal())
				{
					VoronoiVertexInfo temp;
					temp.hDelaunaySimplex=*iter;
					temp.VerticesVectorLocation=vertices.size();
					vInfo.push_back(temp);
					Point p_center=s.center();
					GeometryVector g_center = CGALPoint2GeometryVector(p_center, d);
					vertices.push_back(g_center);
				}

				if(pNeighbors!=nullptr)
					for (int ii = 0; ii <= d; ii++)
					{
						GeometryVector car = CGALPoint2GeometryVector(De.associated_point(De.vertex_of_simplex(*iter, ii)), d);
						GeometryVector rel = c.CartesianCoord2RelativeCoord(car);
						bool found = false;
						size_t n=c.NumParticle();
						c.IterateThroughNeighbors(rel, LengthPrecision*0.01, [&found, &n](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
						{
							n = Sourceparticle;
							found = true;
						}, &found);

						//debug temp
						//std::cout << "found neighbor " << n << '\n';

						if (n < c.NumParticle())
						{
							if (n != i)
								pNeighbors->insert(n);
						}
						else
							std::cerr << "Error in GetVoronoiCell : predicted particle not found!\n";
					}
			}
		}
		//calculate FurthestVertexDist
		FurthestVertexDist=0;
		for(auto iter=vertices.begin(); iter!=vertices.end(); iter++)
		{
			GeometryVector dist=*iter-ci;
			if(FurthestVertexDist<dist.Modulus2())
				FurthestVertexDist=dist.Modulus2();
		}
		FurthestVertexDist=std::sqrt(FurthestVertexDist);

		if(vertices.empty()==false && 2*FurthestVertexDist<radius && pLinkage!=nullptr)
		{
			//calculate linkage
			pLinkage->resize( (d+1)*vertices.size(), -1 );
			std::sort(vInfo.begin(), vInfo.end(), [] (const VoronoiVertexInfo & left, const VoronoiVertexInfo & right) ->bool
			{
				return left.hDelaunaySimplex<right.hDelaunaySimplex;
			});
			for(auto iter=vInfo.begin(); iter!=vInfo.end(); iter++)
			{
				size_t NeighborCount=0;
				for(DimensionType i=0; i<=d; i++)
				{
					Simplex_handle hNeighbor = De.opposite_simplex(iter->hDelaunaySimplex, i);
					auto NewFaceIterator = std::lower_bound(vInfo.begin(), vInfo.end(), hNeighbor, [] (const VoronoiVertexInfo & left, const Simplex_handle & right) ->bool
					{
						return left.hDelaunaySimplex<right;
					});
					if (NewFaceIterator == vInfo.end())
						continue;
					if(NewFaceIterator->hDelaunaySimplex==hNeighbor)
					{
						//debug temp
						//std::cerr<<iter->VerticesVectorLocation<<" \t"<<(d+1)<<" \t"<<NeighborCount<<'\n';
						pLinkage->at( iter->VerticesVectorLocation*(d+1) + NeighborCount ) = NewFaceIterator->VerticesVectorLocation;
						NeighborCount++;
					}
				}
				//for(; NeighborCount<d+1; NeighborCount++)
				//	pLinkage->at( iter->VerticesVectorLocation*(d+1) + NeighborCount ) = -1;
			}
		}
	}
	//debug temp
	//std::cout << "Num iteration = " << NumIteration << '\n';
}

const double Perturbation = 0;
double PolygonVolume_ConvexHull(const std::vector<GeometryVector> & Vertices, DimensionType d, RandomGenerator & gen)
{
	if(Vertices.size()==0)
		return 0;
	ConvexHull De(d);
	GeometryVector center=SameDimensionZeroVector(Vertices[0]);
	for(auto iter=Vertices.begin(); iter!=Vertices.end(); iter++)
	{
		GeometryVector temp(*iter);
		for(DimensionType i=0; i<d; i++)
			temp.x[i]+=Perturbation*temp.x[i]*(gen.RandomDouble()-0.5);
		center.AddFrom(temp);
		cPoint p(d, temp.x, temp.x+d);
		De.insert(p);
	}
	center=center*(1.0/Vertices.size());
	double result=0;
	std::list<ConvexHull::Facet_handle> facets = De.all_facets();
	for(auto iter=facets.begin(); iter!=facets.end(); iter++)
	{
		std::vector<GeometryVector> vecs;
		for(DimensionType i=0; i<d; i++)
		{
			GeometryVector vi=CGALPoint2GeometryVector( (De.point_of_facet(*iter, i)), d);
			vecs.push_back(vi-center);
		}
		result+= ::SimplexVolume(&vecs[0], d);
	}
	return result;
}
double PolygonVolume(const std::vector<GeometryVector> & Vertices, DimensionType d)
{
	if(d==1)
	{
		std::vector<double> coords;
		for(auto iter=Vertices.begin(); iter!=Vertices.end(); iter++)
			coords.push_back( iter->x[0] );
		//std::nth_element(coords.begin(), coords.begin(), coords.end());
		//std::nth_element(coords.begin()+1, coords.end()-1, coords.end());
		//return coords.back()-coords.front();
		return *std::max_element(coords.begin(), coords.end()) - *std::min_element(coords.begin(), coords.end());
	}
	else
	{
		RandomGenerator gen(1);
		return PolygonVolume_ConvexHull(Vertices, d, gen);
	}
}



void PlotVoronoi2D(std::string prefix, const Configuration & List, const std::string & remark)
{
	const Configuration * list= & List;
#ifdef USE_MATHGL

	if(list->GetDimension()==2)
	{
		//calculate the range of the plot
		GeometryVector a0(0, 0);
		GeometryVector a1=a0+list->GetBasisVector(0);
		GeometryVector a2=a0+list->GetBasisVector(1);
		GeometryVector a3=a2+list->GetBasisVector(0);
		double min=0, max=0;
		if(min>a0.x[0])
			min=a0.x[0];
		if(min>a0.x[1])
			min=a0.x[1];
		if(max<a0.x[0])
			max=a0.x[0];
		if(max<a0.x[1])
			max=a0.x[1];
		if(min>a1.x[0])
			min=a1.x[0];
		if(min>a1.x[1])
			min=a1.x[1];
		if(max<a1.x[0])
			max=a1.x[0];
		if(max<a1.x[1])
			max=a1.x[1];
		if(min>a2.x[0])
			min=a2.x[0];
		if(min>a2.x[1])
			min=a2.x[1];
		if(max<a2.x[0])
			max=a2.x[0];
		if(max<a2.x[1])
			max=a2.x[1];
		if(min>a3.x[0])
			min=a3.x[0];
		if(min>a3.x[1])
			min=a3.x[1];
		if(max<a3.x[0])
			max=a3.x[0];
		if(max<a3.x[1])
			max=a3.x[1];

		//prepare the graph
		mglGraph gr;
		double delta=max-min;
		min=min-0.1*delta;
		max=max+0.1*delta;
		gr.SetRanges(min, max, min, max);
		gr.Aspect(1.0, 1.0);
		gr.SetSize(800, 800);
		//gr.SetOrigin(0, 0);
		gr.Puts(mglPoint(0.5*(max+min), min+1.1*(max-min)), remark.c_str());

		//plot crystal BasisVectors
		{
			mglData x(2), y(2);
			x.a[0]=a0.x[0];
			y.a[0]=a0.x[1];
			x.a[1]=a1.x[0];
			y.a[1]=a1.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a1.x[0];
			y.a[0]=a1.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a3.x[0];
			y.a[1]=a3.x[1];
			gr.Plot(x, y, "k-");
			x.a[0]=a2.x[0];
			y.a[0]=a2.x[1];
			x.a[1]=a0.x[0];
			y.a[1]=a0.x[1];
			gr.Plot(x, y, "k-");
		}
		//plot particles
		mglData x(list->NumParticle()), y(list->NumParticle());
		std::vector<double> xx, yy;
		for(int i=0; i<list->NumParticle(); i++)
		{
			GeometryVector now=list->GetCartesianCoordinates(i);
			x.a[i]=now.x[0];
			y.a[i]=now.x[1];

			std::vector<GeometryVector> vertex;
			std::vector<signed long> linkage;
			::GetVoronoiCell(List, i, vertex, &linkage);
			for(size_t i=0; i<linkage.size(); i++)
			{
				if(linkage[i]==-1)
					continue;
				size_t Index1=linkage[i];
				size_t Index2=i/3;
				if(Index1<Index2)
					continue;
				mglData x(2), y(2);
				x.a[0]=vertex[Index1].x[0];
				y.a[0]=vertex[Index1].x[1];
				x.a[1]=vertex[Index2].x[0];
				y.a[1]=vertex[Index2].x[1];
				gr.Plot(x, y, "r-");
			}
			for(auto iter=vertex.begin(); iter!=vertex.end(); iter++)
			{
				xx.push_back(iter->x[0]);
				yy.push_back(iter->x[1]);
			}
		}
		gr.Plot(x, y, "b .");
		mglData xxx(&xx[0], xx.size()), yyy(&yy[0], yy.size());
		gr.Plot(xxx, yyy, "r .");

		prefix += std::string(FigureFormat);
		gr.WriteFrame(prefix.c_str());
		if(AlsoWriteEPS)
		{
			prefix.resize(prefix.size()-4);
			prefix+=std::string(".eps");
			gr.WriteEPS(prefix.c_str());
		}
	}
	else
	{
		std::cerr<<"error in PlotVoronoi2D : unsupported dimension\n";
	}
#endif
}
void DisplayConfiguration_3D_WithVoronoi(const Configuration & list, double Radius, long Num)
{
	assert(list.GetDimension()==3);
	long & TargetI=Num;
	//debug temp
	//std::cin>>TargetI;
	double typicalLength=std::pow(list.PeriodicVolume()/list.NumParticle(), 1.0/3.0);
	for(size_t i=0; i<list.NumParticle(); i++)
	{
		//Configuration::particle * pa=list.GetParticle(i);
		list.IterateThroughNeighbors(i, typicalLength, [&typicalLength, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
		{
			if(SourceAtom!=i)
			{
				double l=std::sqrt(shift.Modulus2());
				if(l<typicalLength)
					typicalLength=l;
			}
		});
	}

	std::vector<sphere> spheres;
	std::vector<line> lines;

	std::vector<GeometryVector> a(8, GeometryVector(0.0, 0.0, 0.0));
	a[1]=a[0]+list.GetBasisVector(0);
	a[2]=a[0]+list.GetBasisVector(1);
	a[3]=a[0]+list.GetBasisVector(2);
	a[4]=a[1]+list.GetBasisVector(1);
	a[5]=a[1]+list.GetBasisVector(2);
	a[6]=a[2]+list.GetBasisVector(2);
	a[7]=a[4]+list.GetBasisVector(2);
	GeometryVector shift=a[7]*0.5;
	double one_over_scale=1.0/std::sqrt(shift.Modulus2());
	for(int i=0; i<8; i++)
	{
		a[i].MinusFrom(shift);
		a[i].MultiplyFrom(one_over_scale);
	}
	int point1[12]={0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6};
	int point2[12]={1, 2, 3, 4, 5, 4, 6, 5, 6, 7, 7, 7};
	for(int i=0; i<12; i++)
	{
		line temp;
		temp.x1=a[point1[i]].x[0];
		temp.y1=a[point1[i]].x[1];
		temp.z1=a[point1[i]].x[2];
		temp.x2=a[point2[i]].x[0];
		temp.y2=a[point2[i]].x[1];
		temp.z2=a[point2[i]].x[2];
		temp.blue=0.0;
		temp.red=0.0; 
		temp.green=0.0;
		temp.transparency=1.0;
		temp.width = 0.06*typicalLength*one_over_scale;
		lines.push_back(temp);
	}

	for(size_t i=0; i<list.NumParticle(); i++)
	{
		//auto pA=list.GetParticle(i);
		GeometryVector loc=(list.GetCartesianCoordinates(i)-shift)*one_over_scale;
		sphere temp;
		temp.x=loc.x[0];
		temp.y=loc.x[1];
		temp.z=loc.x[2];

		char chr1=list.GetCharacteristics(i).name[0];
		if(chr1>='A' && chr1<='Z')
		{
			temp.red=((chr1-'A')/6)*0.1+0.3;
			temp.green=((chr1-'A')%6)*0.1+0.3;
			char chr2=list.GetCharacteristics(i).name[1];
			if(chr2>='a' && chr2<='z')
				temp.blue=(chr2-'a')*0.02+0.4;
			else
				temp.blue=0.35;
		}
		else
		{
			temp.red=0.25;
			temp.green=0.25;
			temp.blue=0.25;
		}
		temp.transparency=0.97;

		if(Radius==0.0)
			temp.radius=0.3*typicalLength*one_over_scale;
		else
			temp.radius=Radius*one_over_scale;

		spheres.push_back(temp);

		if(i==TargetI || TargetI==-1)
		{
			std::vector<GeometryVector> vertex;
			std::vector<signed long> linkage;
			::GetVoronoiCell(list, i, vertex, &linkage);
			for(size_t i=0; i<linkage.size(); i++)
			{
				if(linkage[i]==-1)
					continue;
				size_t Index1=linkage[i];
				size_t Index2=i/4;
				if(Index1<Index2)
					continue;
				GeometryVector end1=(vertex[Index1]-shift)*one_over_scale;
				GeometryVector end2=(vertex[Index2]-shift)*one_over_scale;
				line temp;
				temp.x1=end1.x[0];
				temp.x2=end2.x[0];
				temp.y1=end1.x[1];
				temp.y2=end2.x[1];
				temp.z1=end1.x[2];
				temp.z2=end2.x[2];
				temp.blue=1.0;
				temp.red=0.0; 
				temp.green=0.0;
				temp.transparency=1.0;
				temp.width = 0.03*typicalLength*one_over_scale;
				lines.push_back(temp);
			}
			for(size_t i=0; i<vertex.size(); i++)
			{
				GeometryVector loc=(vertex[i]-shift)*one_over_scale;
				sphere temp;
				temp.x=loc.x[0];
				temp.y=loc.x[1];
				temp.z=loc.x[2];

				temp.blue=1.0;
				temp.red=0.0; 
				temp.green=0.0;
				temp.transparency=0.97;

				if(Radius==0.0)
					temp.radius=0.1*typicalLength*one_over_scale;
				else
					temp.radius=Radius*one_over_scale*0.33;

				spheres.push_back(temp);
			}
		}

	}

	DisplaySpheres(spheres, lines);
}
#include <omp.h>
void VoronoiVolumeDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleSize, double ResolutionPreference)
{
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"Generate VoronoiVolumeDistrubution";
	progress_display pd(NumConfigs);

	std::vector<double> distances;

	DimensionType dim;

	size_t NumConfigsSampled=0;

	omp_lock_t lock;
	omp_init_lock(&lock);

	auto GetVoronoiSizeFunction = [] (const Configuration & c, size_t i, DimensionType dim) ->double
	{
		std::vector<GeometryVector> vertices;
		::GetVoronoiCell(c, i, vertices);
		return PolygonVolume(vertices, dim);
	};
	for(size_t j=0; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		if(j==0)
			dim=Config.GetDimension();
		size_t NumParticle=Config.NumParticle();
		double TypicalLength=std::pow(Config.PeriodicVolume()/Config.NumParticle(), 1.0/Config.GetDimension())*1.5;

		long end=Config.NumParticle();
#pragma omp parallel 
		{
			Configuration ThreadConfig(Config);
#pragma omp for schedule(dynamic)
			for(long j=0; j<end; j++)
			{
				if(Verbosity>4)
					std::cout<<"Doing Particle"<<j<<'\n';
				double result=GetVoronoiSizeFunction(ThreadConfig, j, dim);
				omp_set_lock(&lock);
				distances.push_back(result);
				omp_unset_lock(&lock);
			}
		}
		NumConfigsSampled++;

		if(distances.size()>SampleSize)
			break;
	}
	omp_destroy_lock(&lock);

	//if(Verbosity>2)
	//	std::cout<<"generating bins...";
	//prepare for the iteration
	Result.clear();
	//iteration
	HistogramGenerator HGen;
	HGen.GenerateBins(distances, static_cast<double>(NumConfigs)/NumConfigsSampled, ResolutionPreference);

	//Re-do the bin counts
	//if(Verbosity>2)
	//	std::cout<<"counting pair distances...";

	for(size_t j=NumConfigsSampled; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		size_t NumParticle=Config.NumParticle();

		long end=Config.NumParticle();
#pragma omp parallel 
		{
			Configuration ThreadConfig(Config);
#pragma omp for schedule(dynamic)
			for(long j=0; j<end; j++)
			{
				if(Verbosity>4)
					std::cout<<"Doing Particle"<<j<<'\n';
				HGen.Report(GetVoronoiSizeFunction(ThreadConfig, j, dim));
			}
		}
	}

	size_t TotalCount=0;
	for(auto iter=HGen.bins.begin(); iter!=HGen.bins.end(); iter++)
		TotalCount+=iter->Count;

	for(auto iter=HGen.bins.begin(); iter!=HGen.bins.end(); iter++)
	{
		GeometryVector temp(4);
		if(iter!=HGen.bins.end()-1)
			temp.x[0]=((iter[0].Start)+(iter[1].Start))/2.0;//r 
		else
			continue;

		if(iter==HGen.bins.begin())
		{
			if(iter->Start==0.0 && (iter[2].Start-iter[1].Start)>iter[1].Start && iter->Count==0)
				continue;//there shouldn't be a bin at here
		}
		temp.x[2]=(temp.x[0]-iter[0].Start);//\delta r

		temp.x[1]=static_cast<double>(iter->Count)/TotalCount/(2*temp.x[2]);//p(r)

		if(iter->Count!=0)
			temp.x[3]=temp.x[1]/std::sqrt(static_cast<double>(iter->Count));//\delta p(r)
		else
			temp.x[3]=0.0;
		Result.push_back(temp);
	}

	if(Verbosity>2)
		std::cout<<"done!\n";
}
void VoronoiNumSidesDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig, std::vector<GeometryVector> & Result)
{
	//get the pair distances list
	if (Verbosity>2)
		std::cout << "Generate VoronoiNumSidesDistrubution";
	std::map<size_t, size_t> sides;
	progress_display pd(NumConfig);
	omp_lock_t pdlock;
	omp_init_lock(&pdlock);
	omp_lock_t sideslock;
	omp_init_lock(&sideslock);
	size_t SumNumParticle = 0;
#pragma omp parallel for
	for (signed long i = 0; i < NumConfig; i++)
	{
		Configuration c = GetConfigsFunction(i);
		size_t np = c.NumParticle();
#pragma omp atomic
		SumNumParticle += np;
		for (size_t j = 0; j < np; j++)
		{
			std::vector<GeometryVector> vertices;
			std::vector<signed long> linkage;
			GetVoronoiCell(c, j, vertices, &linkage);
			size_t nside = 0;
			for (auto iter = linkage.begin(); iter != linkage.end(); ++iter)
				if (*iter != -1)
					nside++;
			assert(nside % 2 == 0);
			nside /= 2;
			omp_set_lock(&sideslock);
			std::map<size_t, size_t>::iterator iter = sides.find(nside);
			if (iter == sides.end())
				sides.insert(std::make_pair(nside, (size_t)(1)));
			else
				iter->second++;
			omp_unset_lock(&sideslock);
		}
		omp_set_lock(&pdlock);
		pd++;
		omp_unset_lock(&pdlock);
	}
	omp_destroy_lock(&sideslock);
	omp_destroy_lock(&pdlock);
	for (auto iter = sides.begin(); iter != sides.end(); iter++)
		Result.push_back(GeometryVector(iter->first, (double)(iter->second) / SumNumParticle));

	if (Verbosity>2)
		std::cout << "done!\n";
}


//C(r)=<(V(x)-<V>)(V(x+r)-<V>)>_x
void VoronoiVolumeCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, double MaxDistance, double RPrecision)
{
	//get the pair distances list
	if(Verbosity>2)
		std::cout<<"Generate VoronoiVolumeCorrelation";
	progress_display pd(NumConfigs);

	std::vector<double> distances;

	DimensionType dim;

	auto GetVoronoiSizeFunction = [] (const Configuration & c, size_t i, DimensionType dim) ->double
	{
		std::vector<GeometryVector> vertices;
		::GetVoronoiCell(c, i, vertices);
		return PolygonVolume(vertices, dim);
	};
	//double SumVolume=0;
	//size_t SumNumParticle=0;

	//this function calculates <(V(x)-<V>)(V(x+r)-<V>)>_x
	//this struct is a single data of (V(x)-<V>)(V(x+r)-<V>)
	struct DataType
	{
		double r;
		double vv;
	};
	std::vector<DataType> vData;

	for(size_t j=0; j<NumConfigs; j++)
	{
		pd++;
		//if(Verbosity>3 || (Verbosity>2&&j%100==0) )
		//	std::cout<<j<<"/"<<NumConfigs<<"configurations processed\n";
		Configuration Config = GetConfigsFunction(j);
		if(j==0)
			dim=Config.GetDimension();
		size_t NumParticle=Config.NumParticle();
		//SumNumParticle+=NumParticle;
		double AverageVolume=Config.PeriodicVolume()/Config.NumParticle();
		std::vector<double> VVolume(NumParticle, 0.0);

		long end=NumParticle;
#pragma omp parallel 
		{
			Configuration ThreadConfig(Config);
#pragma omp for schedule(dynamic)
			for(long j=0; j<end; j++)
			{
				if(Verbosity>4)
					std::cout<<"Doing Particle"<<j<<'\n';
				double result=GetVoronoiSizeFunction(ThreadConfig, j, dim);
				VVolume[j]=result-AverageVolume;
//#pragma omp atomic
//				SumVolume+=result;
			}
		}

		for(long j=0; j<end; j++)
		{
			auto func=[&] (const GeometryVector &shift, const GeometryVector &LatticeShift, const signed long *PeriodicShift, const size_t SourceAtom) ->void
			{
				DataType temp;
				if(SourceAtom==j && shift.Modulus2()<1e-10)
					temp.r=0.0;
				else
					temp.r=std::sqrt(shift.Modulus2());
				temp.vv=VVolume[j]*VVolume[SourceAtom];
				vData.push_back(temp);
			};
			GeometryVector Rel=Config.GetRelativeCoordinates(j);
			Config.IterateThroughNeighbors(Rel, MaxDistance, func);
		}
	}

//	double AverageVolume=SumVolume/SumNumParticle;
	auto CompareFunc=[] (const DataType & left, const DataType & right) ->bool
	{
		return left.r<right.r;
	};
	std::sort(vData.begin(), vData.end(), CompareFunc);

	size_t SameRCount=0;
	double Sum=0.0;
	double Sum2=0.0;
	double CurrentR=0.0;
	double SumR2=0.0;
	//if(KPrecision==0.0)
	//	KPrecision = LengthPrecision/tempList.GetBasisVector(0).Modulus2();
	for(auto iter=vData.begin(); iter!=vData.end(); iter++)
	{
		double c=iter->vv;
		double c2=c*c;
		double r=iter->r;
		if( (r-CurrentR < RPrecision&& CurrentR>0.0) || (r==0.0 && CurrentR==0.0) )
		{
			SameRCount++;
			Sum+=c;
			Sum2+=c2;
			SumR2+=r*r;
		}
		else
		{
			if(SameRCount>0)
			{
				GeometryVector temp(4);
				temp.x[0]=std::sqrt(SumR2/SameRCount);
				temp.x[1]=Sum/SameRCount;
				if(CurrentR>0.0)
					temp.x[2]=RPrecision;
				else
					temp.x[2]=0;
				temp.x[3]=std::sqrt(Sum2/(SameRCount)-temp.x[1]*temp.x[1])/(SameRCount);
				Result.push_back(temp);
			}

			Sum=c;
			Sum2=c*c;
			SameRCount=1;
			CurrentR=r;
			SumR2=r*r;
		}
	}
	if(SameRCount>0)
	{
		GeometryVector temp(4);
		temp.x[0]=std::sqrt(SumR2/SameRCount);
		temp.x[1]=Sum/SameRCount;
		temp.x[2]=RPrecision;
		temp.x[3]=std::sqrt(Sum2/(SameRCount)-temp.x[1]*temp.x[1])/(SameRCount);
		Result.push_back(temp);
	}

	if(Verbosity>2)
		std::cout<<"done!\n";
}

#else

void GetVoronoiCell(const PeriodicCellList<Empty> & c, size_t i, std::vector<GeometryVector> & vertices, std::vector<signed long> * pLinkage, std::set<size_t> * pNeighbors)
{
	std::cerr<<"Error : CGAL not enabled!\n";
	return;
}
double PolygonVolume(std::vector<GeometryVector> Vertices, DimensionType d)
{
	std::cerr<<"Error : CGAL not enabled!\n";
	return 0;
}
void PlotVoronoi2D(std::string prefix, const Configuration & List, const std::string & remark)
{
	std::cerr<<"Error : CGAL not enabled!\n";
	return;
}

//Num=-1 : display voronoi of all particles Num>=0 : display voronoi of a specific particle
void DisplayConfiguration_3D_WithVoronoi(const Configuration & list, double Radius, long Num)
{
	std::cerr<<"Error : CGAL not enabled!\n";
	return;
}

void VoronoiVolumeDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, size_t SampleSize, double ResolutionPreference)
{
	std::cerr<<"Error : CGAL not enabled!\n";
	return;
}
void VoronoiVolumeCorrelation(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfigs, std::vector<GeometryVector> & Result, double MaxDistance, double RPrecision)
{
	std::cerr<<"Error : CGAL not enabled!\n";
	return;
}

void VoronoiNumSidesDistrubution(std::function<const Configuration(size_t i)> GetConfigsFunction, size_t NumConfig, std::vector<GeometryVector> & Result)
{
	std::cerr<<"Error : CGAL not enabled!\n";
	return;
}
#endif