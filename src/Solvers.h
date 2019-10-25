#ifndef SOLVERS_INCLUDED
#define SOLVERS_INCLUDED


#include "LatticeSumSolver.h"
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

class FromStructure : public LatticeSumSolver
{
public:
	Configuration Structure;
	FromStructure(const Configuration & source):Structure(source)
	{
	}
	virtual const char * Tag(void)
	{
		return "TempSolver";
	}
	virtual Configuration GetStructure(void)
	{
		return this->Structure;
	}
};

class LoadStructure : public LatticeSumSolver
{
public:
	std::string filename;
	Configuration * pStructure;
	LoadStructure(const char * FileName):filename(FileName)
	{
		std::fstream ifile(this->filename.c_str(), std::fstream::in);

		if(!ifile.good())
		{
			std::cerr<<"error in LoadStructure::GetStructure : file"<<filename<<"is invalid!\n";
			exit(1);
		}
		this->pStructure= new Configuration(::ReadPos(ifile));
	}
	LoadStructure(const LoadStructure & src) : filename(src.filename), pStructure(new Configuration(*src.pStructure))
	{
	}
	virtual const char * Tag(void)
	{
		return filename.c_str();
	}
	virtual Configuration GetStructure(void)
	{
		return *this->pStructure;
	}
	virtual ~LoadStructure()
	{
		delete this->pStructure;
	}
};

class CloneSolver : public LatticeSumSolver
{
private:
	std::string tag;
	Configuration stru;
public:
	CloneSolver(LatticeSumSolver & source) : stru(source.GetStructure())
	{
		this->tag.assign(source.Tag());
		this->Terms.assign(source.Terms.begin(), source.Terms.end());
	}
	virtual const char * Tag(void)
	{
		return tag.c_str();
	}
	virtual Configuration GetStructure(void)
	{
		return this->stru;
	}
};

class twoDSquare : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "twoDSquare";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0));
		base.push_back(GeometryVector(0, 1));
		Configuration result(2, &base[0], 1.0);

		result.Insert("C", GeometryVector(0, 0));

		return result;
	}
};
class twoDTriangle : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "twoDTriangle";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0));
		base.push_back(GeometryVector(0.5, std::sqrt((double)(3))*0.5));
		Configuration result(2, &base[0], 1.0);

		result.Insert("C", GeometryVector(0, 0));


		result.Rescale(1+ 1e-11);

		return result;
	}
};
class twoDKagome : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "twoDKagome";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(2, 0));
		base.push_back(GeometryVector(1, std::sqrt((double)(3))));
		Configuration result(2, &base[0], 1.0);

		result.Insert("C", GeometryVector(0.5, 0));
		result.Insert("C", GeometryVector(0, 0.5));
		result.Insert("C", GeometryVector(0.5, 0.5));

		return result;
	}
};
class twoDGraphene : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "twoDGraphene";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0));
		base.push_back(GeometryVector(-0.5, std::sqrt((double)(3))*0.5));
		Configuration result(2, &base[0], 1.0);

		result.Insert("C", GeometryVector(0, 0));
		result.Insert("C", GeometryVector(1.0/3, 2.0/3));

		return result;
	}
};

class Graphite : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "graphite";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0, 0));
		base.push_back(GeometryVector(0.5, std::sqrt((double)(3))*0.5, 0));
		base.push_back(GeometryVector(0,0,2.726));
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0, 0, 0));
		result.Insert("C", GeometryVector(0, 0, 0.5));
		result.Insert("C", GeometryVector(1.0/3, 1.0/3, 0));
		result.Insert("C", GeometryVector(2.0/3, 2.0/3, 0.5));

		return result;
	}
};

class SimpleCubic : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "simple cubic";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0, 0));
		base.push_back(GeometryVector(0,1,0));
		base.push_back(GeometryVector(0,0,1));
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0, 0, 0));

		return result;
	}
};

class FCC : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "fcc";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(0.5,0.5,0));
		base.push_back(GeometryVector(0,0.5,0.5));
		base.push_back(GeometryVector(0.5,0,0.5));
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0, 0, 0));

		return result;
	}
};
class HCP : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "HCP";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0, 0));
		base.push_back(GeometryVector(0.5, std::sqrt((double)(3))*0.5, 0));
		base.push_back(GeometryVector(0 , 0, std::sqrt((double)(8)/3)));
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0, 0, 0));
		result.Insert("C", GeometryVector((double)(1)/3, (double)(1)/3, 0.5));

		return result;
	}
};
class Diamond : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "diamond";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(0.5,0.5,0));
		base.push_back(GeometryVector(0,0.5,0.5));
		base.push_back(GeometryVector(0.5,0,0.5));
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0, 0, 0));
		result.Insert("C", GeometryVector(0.25,0.25,0.25));

		return result;
	}
};
class Kagome : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "kagome";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(0.5,0.5,0));
		base.push_back(GeometryVector(0,0.5,0.5));
		base.push_back(GeometryVector(0.5,0,0.5));
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0.125,0.125,0.125));
		result.Insert("C", GeometryVector(0.625,0.125,0.125));
		result.Insert("C", GeometryVector(0.125,0.625,0.125));
		result.Insert("C", GeometryVector(0.125,0.125,0.625));

		return result;
	}

};

class BCC : public LatticeSumSolver
{
public:
	virtual const char * Tag(void)
	{
		return "bcc";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0, 0));
		base.push_back(GeometryVector(0, 1, 0));
		base.push_back(GeometryVector(0, 0, 1));
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0,0,0));
		result.Insert("C", GeometryVector(0.5,0.5,0.5));

		return result;
	}
};

class GeneralizedGraphite : public LatticeSumSolver
{
public:
	GeometryVector vec2;
	GeneralizedGraphite(const GeometryVector & BaseVector2):vec2(BaseVector2)
	{
	}
	virtual const char * Tag(void)
	{
		return "GeneralizedGraphite";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0, 0));
		base.push_back(GeometryVector(0.5, std::sqrt((double)(3))*0.5, 0));
		base.push_back(this->vec2);
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0,0,0));
		result.Insert("C", GeometryVector(1.0/3,1.0/3,0));

		return result;
	}

};
class GeneralizedHexagonal : public LatticeSumSolver
{
public:
	GeometryVector vec2;
	GeneralizedHexagonal(const GeometryVector & BaseVector2):vec2(BaseVector2)
	{
	}
	virtual const char * Tag(void)
	{
		return "GeneralizedHexagonal";
	}
	virtual Configuration GetStructure(void)
	{
		std::vector<GeometryVector> base;
		base.push_back(GeometryVector(1, 0, 0));
		base.push_back(GeometryVector(0.5, std::sqrt((double)(3))*0.5, 0));
		base.push_back(this->vec2);
		Configuration result(3, &base[0], 1.0);

		result.Insert("C", GeometryVector(0,0,0));

		return result;
	}
};


#endif