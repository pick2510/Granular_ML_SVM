#include "PeriodicCellList.h"
#include <exception>



void Output(std::string prefix, const Configuration & List)
{
	prefix+=std::string(".pos");
	std::fstream out(prefix.c_str(), std::fstream::out);
	out.precision(17);

	Output(out, List);
}

void Output(std::ostream & out, const Configuration & List)
{
	const Configuration * list = & List;
	//determine Rc of output
	//GeometryVector temp(list->GetDimension());
	//for( ::DimensionType i=0; i<list->GetDimension(); i++)
	//{
	//	GeometryVector t=list->GetBasisVector(i);
	//	if(temp.Dot(t)>0)
	//		temp.AddFrom(t);
	//	else
	//		temp.MinusFrom(t);
	//}
	//double Rc=2*std::sqrt(temp.Modulus2())/std::pow(list->NumParticle(), 1.0/list->GetDimension());
	double Rc=std::pow(100*list->PeriodicVolume()/list->NumParticle(), 1.0/list->GetDimension());
	//////////////////////////////////////
	GeometryVector origin(list->GetDimension());
	size_t count=0;
	list->IterateThroughNeighbors(origin, Rc, [&count](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void {count++;} );
	out<<count<<'\n';
	//the second line is a comment
	out<<'\n';
	//.xyz part, used by visualization sofwares
	list->IterateThroughNeighbors(origin, Rc, [&out, &list](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom) ->void 
	{
		out<<list->GetCharacteristics(SourceAtom).name<<" \t"; 
		shift.OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	} );


	out<<"\n\n****************************************\n\n";

	out<<"Primitive vectors\n";
	for(::DimensionType i=0; i<list->GetDimension(); i++)
	{
		out<<"a(";
		out<<i;
		out<<") = \t";
		list->GetBasisVector(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}

	out<<"\n\nVolume =  "<<list->PeriodicVolume()<<"\n\n";
	out<<"Reciprocal vectors\n";
	for(::DimensionType i=0; i<list->GetDimension(); i++)
	{
		out<<"b(";
		out<<i;
		out<<") = \t";
		list->GetReciprocalBasisVector(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}
	out<<"\n\nBasis Vectors:\nAtom    Lattice Coordinates                Cartesian Coordinates\n\n";
	for(size_t i=0; i<list->NumParticle(); i++)
	{
		//Configuration::particle * now=list->GetParticle(i);
		out<<list->GetCharacteristics(i).name<<" \t";
		list->GetRelativeCoordinates(i).OutputCoordinate(out, list->GetDimension());
		list->GetCartesianCoordinates(i).OutputCoordinate(out, list->GetDimension());
		out<<'\n';
	}
}
Configuration ReadPos(std::string prefix)
{
	std::string name=prefix+std::string(".pos");
	std::fstream file(name.c_str(), std::fstream::in);
	return ReadPos(file);
}
Configuration ReadPos(std::istream & ifile)
{
	const std::string flag1("Primitive");
	const std::string flag2("Basis");

	std::string current;

	size_t NumParticle=0;

	::DimensionType dimension=0;
	do
	{
		std::getline(ifile, current);
		if(ifile.eof() || ifile.fail())
		{
			return Configuration(0, nullptr, 0.0, false);
		}
	}
	while(current.find(flag1)!=0);
	for(;;)
	{
		std::getline(ifile, current);
		if(current.size()<3)
			break;
		if(current[2]>='0' && current[2]<='9')
			dimension++;
		else
			break;
	}

	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag2)!=0);
	std::getline(ifile, current);
	std::getline(ifile, current);
	while(ifile.eof()==false)
	{
		char junk[1000];
		ifile.getline(junk, 1000);

		if(ifile.eof()==false)
			NumParticle++;
	}
	ifile.clear();
	ifile.seekg(0, std::ios::beg);


	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag1)!=0);

	std::vector<GeometryVector> base;
	for(int i=0; i<dimension; i++)
	{
		char a;
		do
		{
			a=ifile.get();
		}
		while(a!='=');
		base.push_back(GeometryVector(dimension));
		base.back().InputCoordinate(ifile, dimension);
	}
	double volume= ::Volume(&base[0], dimension);
	double cellsize = std::pow(volume/(NumParticle/ParticlePerCell), static_cast<double>(1)/dimension);
	Configuration result(dimension, &base[0], cellsize);

	do
	{
		std::getline(ifile, current);
	}
	while(current.find(flag2)!=0);
	std::getline(ifile, current);
	std::getline(ifile, current);
	while(ifile.eof()==false)
	{
		char name[128];
		ifile >> name;

		GeometryVector loc(dimension);
		loc.InputCoordinate(ifile, dimension);

		GeometryVector junk(dimension);
		junk.InputCoordinate(ifile, dimension);

		if(ifile.eof()==false)
			result.Insert(name, loc);
	}
	return result;
}
Configuration ReadHoomdXml(std::string prefix)
{
	prefix+=std::string(".xml");
	std::fstream ifile(prefix.c_str());
	return ReadHoomdXml(ifile);
}
Configuration ReadHoomdXml(std::istream & ifile)
{
	std::string temp;
	std::string flag1("<configuration");
	std::string flag2("<type");

	DimensionType dim;
	size_t NumParticle;
	double volume=1;

	std::vector<GeometryVector> basis;
	std::vector<double> coords;

	do
	{
		std::getline(ifile, temp);
	}
	while(temp.find(flag1)!=0);

	{
		std::stringstream stream;
		stream<<temp;

		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>NumParticle;//this is actually timestep. junk it
		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>dim;
		//read a number
		for(;;)
		{
			char c=stream.peek();
			if(c>='0' && c<='9')
				break;
			if(c=='.')
				break;
			if(stream.eof())
				break;
			stream>>c;
		}
		stream>>NumParticle;
	}

	std::getline(ifile, temp);
	{
		std::stringstream stream;
		stream<<temp;
		double length;

		for(DimensionType d=0; d<dim; d++)
		{
			//read a number
			for(;;)
			{
				char c=stream.peek();
				if(c>='0' && c<='9')
					break;
				if(c=='.')
					break;
				if(stream.eof())
					break;
				stream>>c;
			}
			stream>>length;
			GeometryVector tempVec(dim);
			tempVec.x[d]=length;
			basis.push_back(tempVec);
			volume=volume*length;
		}
	}

	double cellsize = std::pow(volume/(NumParticle/ParticlePerCell), static_cast<double>(1)/dim);
	Configuration result(dim, &basis[0], cellsize);
	std::getline(ifile, temp);
	for(size_t i=0; i<dim*NumParticle; i++)
	{
		double tempD;
		ifile>>tempD;
		coords.push_back(tempD);
	}


	do
	{
		std::getline(ifile, temp);
	}
	while(temp.find(flag2)!=0);
	for(size_t i=0; i<NumParticle; i++)
	{
		std::getline(ifile, temp);
		assert(temp.length()<3);
		GeometryVector rel(dim);
		for(DimensionType j=0; j<dim; j++)
			rel.x[j]=coords[i*dim+j]/basis[j].x[j];

		result.Insert(temp.c_str(), rel);
	}

	return result;
}
Configuration ReadCoordinate(std::istream & file)
{
	return ReadCoordinate(file, 1.0);
}
Configuration ReadCoordinate(std::istream & file, double SideLength)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;

	//determine dimension
	std::getline(file, temp);
	std::stringstream str;
	str<<temp;
	for(;;)
	{
		double dtemp=-1.0;
		str>>dtemp;
		if(dtemp!=-1.0)
			dim++;
		else
			break;
	}

	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
			file>>rel.x[j];
		result.Insert("A", rel);
	}

	return result;
}
Configuration ReadCartesianCoordinate(std::istream & file)
{
	return ReadCartesianCoordinate(file, 1.0);
}
Configuration ReadCartesianCoordinate(std::istream & file, double SideLength, DimensionType dimension)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;

	if (dimension == 0)
	{
		//determine dimension
		std::getline(file, temp);
		std::stringstream str;
		str << temp;
		for (;;)
		{
			double dtemp = -1.0;
			str >> dtemp;
			if (dtemp != -1.0)
				dim++;
			else
				break;
		}
	}
	else
		dim = dimension;

	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for (;;)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
			file>>rel.x[j];
		if (file.good())
		{
			rel = result.CartesianCoord2RelativeCoord(rel);
			result.Insert("A", rel);
		}
		else
			break;
	}

	return result;
}
Configuration ReadGro(std::istream & file)
{
	GeometryVector basis[ ::MaxDimension];
	DimensionType dim=3;
	size_t NumP=0;
	std::string temp;

	//determine NumP
	std::getline(file, temp);
	file>>NumP;
	for(int i=0; i<NumP; i++)
	{
		std::getline(file, temp);
	}
	std::getline(file, temp);
	//determine box size
	double Volume=1.0;
	for(int i=0; i<3; i++)
	{
		double l;
		file>>l;
		basis[i]=GeometryVector(3);
		basis[i].x[i]=l;
		Volume*=l;
	}
	file.clear();
	file.seekg(0, std::ios::beg);
	std::getline(file, temp);
	std::getline(file, temp);

	double CellSize=Volume/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		//Only read atoms named 'O'!
		file>>temp;
		file>>temp;
		if(temp.at(0)=='O')
		{
			file>>temp;
			GeometryVector rel(dim);
			for(size_t j=0; j<dim; j++)
				file>>rel.x[j];
			rel=result.CartesianCoord2RelativeCoord(rel);
			result.Insert("A", rel);
		}
		else
			std::getline(file, temp);
	}

	return result;
}
Configuration ReadStealthOutput(std::istream & file, double SideLength)
{
	DimensionType dim=0;
	size_t NumP=0;
	std::string temp;
	std::vector<double> OneOverSize;//the side length of the input file

	//determine dimension
	std::getline(file, temp);
	{
		std::stringstream str;
		str<<temp;
		str>>temp;
		str>>dim;
	}

	//determine size
	std::getline(file, temp);
	{
		std::stringstream str;
		str<<temp;
		str>>temp;
		for(DimensionType i=0; i<dim; i++)
		{
			double tt;
			str>>tt;
			OneOverSize.push_back(1.0/tt);
		}
	}

	std::getline(file, temp);
	//determine NumP
	while(file.eof()==false)
	{
		std::getline(file, temp);
		NumP++;
	}
	file.clear();
	file.seekg(0, std::ios::beg);
	std::getline(file, temp);
	std::getline(file, temp);

	std::vector<GeometryVector> basis;
	for(DimensionType i=0; i<dim; i++)
	{
		basis.push_back(GeometryVector(dim));
		basis.back().x[i]=SideLength;
	}
	double CellSize=SideLength/std::pow(static_cast<double>(NumP), 1.0/static_cast<double>(dim));
	Configuration result(dim, &basis[0], CellSize);

	for(size_t i=0; i<NumP; i++)
	{
		GeometryVector rel(dim);
		for(size_t j=0; j<dim; j++)
		{
			file>>rel.x[j];
			rel.x[j]*=OneOverSize[j];
		}
		result.Insert("A", rel);
	}

	return result;
}
Configuration ReadStealthOutput(std::istream & file)
{
	return ReadStealthOutput(file, 1.0);
}

SpherePacking ReadTJOutput(std::istream & inFile)
{
	std::string input;
	std::getline(inFile,input);
	std::stringstream ss1(input);
	int dim, N;
	ss1 >> dim;
	std::string buffer;
	ss1 >> buffer; //"HS"
	std::string packType;
	ss1 >> packType; //"poly" or "mono"
	getline(inFile,input); //Line 2 is junk.
	inFile >> N;

	double radMono;
	if (packType == "mono") 
	{//This line is the diameter value for monodisperse packings
		inFile >> radMono;
		radMono /= 2.0; //Convert from diameter to radius
	}
	else if (packType != "poly")
	{
		std::cerr << "Error in ReadTJOutput : Unsupported packType!\n";
		return SpherePacking();
	}

	// lattice vectors next
	GeometryVector basis[ ::MaxDimension];
	for (int i=0; i<dim; i++) 
	{
		basis[i]=GeometryVector(dim);
		for (int j=0; j<dim; j++) 
		{
			double val;
			inFile >> val;
			basis[i].x[j]=val;
		}
	}
	SpherePacking packing(dim, basis, ::MaxDistance);

	// then more trash
	getline(inFile,input);
	getline(inFile,input);

	// then coordinates and radii
	for (int i=0; i<N; i++) 
	{
		double radii;
		GeometryVector car(dim);
		for (int d=0; d<dim; d++) 
		{
			inFile >> car.x[d];
		}
		if (packType == "poly") 
		{
			inFile >> radii;
		}
		else 
		{
			radii = radMono;
		}
		packing.Insert(radii, packing.CartesianCoord2RelativeCoord(car));
	}
	packing.SetCellSize(packing.GetMaxRadius()*1.05);

	return packing; //Got a packing!

}
SpherePacking ReadTJOutput(const std::string & prefix)
{
	std::string name = prefix + std::string(".dat");
	std::fstream ifile(name.c_str(), std::fstream::in);
	return ReadTJOutput(ifile);
}
void WriteTJOutput(std::fstream & ofile, const SpherePacking & pk)
{
	ofile.precision(17);
	DimensionType dim = pk.GetDimension();
	size_t num = pk.NumParticle();
	if (dim == 0 || num == 0)
	{
		std::cerr << "Error in WriteTJOutput : invalid packing!\n";
		return;
	}
	//determine mono or poly
	bool mono = true;
	for(size_t i=1; i<num; i++)
		if (pk.GetCharacteristics(i) != pk.GetCharacteristics(0))
		{
			mono = false;
			break;
		}

	ofile << dim << "	hs	";
	if (mono)
		ofile << "mono";
	else
		ofile << "poly";
	ofile << "\n" << num << " 1\n" << num << '\n';
	if (mono)
		ofile << 2 * pk.GetCharacteristics(0) << '\n';

	for (DimensionType i = 0; i < dim; i++)
	{
		for (DimensionType j = 0; j < dim; j++)
			ofile << pk.GetBasisVector(i).x[j] << " \t";
		ofile << '\n';
	}
	ofile << "T T T";
	for (size_t n = 0; n < num; n++)
	{
		ofile << "\n";
		for (DimensionType j = 0; j < dim; j++)
			ofile << pk.GetCartesianCoordinates(n).x[j] << " \t";
		if (!mono)
			ofile << pk.GetCharacteristics(n);
	}
}
void WriteTJOutput(const std::string & prefix, const SpherePacking & pk)
{
	std::string name = prefix + std::string(".dat");
	std::fstream ifile(name.c_str(), std::fstream::out);
	return WriteTJOutput(ifile, pk);
}


Configuration MultiplicateStructure(const Configuration & src, size_t NumberInEachSide)
{
	double dNumberInEachSide=static_cast<double>(NumberInEachSide);
	DimensionType dim=src.GetDimension();
	std::vector<GeometryVector> bases;
	for(DimensionType i=0; i<dim; i++)
		bases.push_back(src.GetBasisVector(i)*dNumberInEachSide);
	Configuration result(src);
	result.ChangeBasisVector(&bases[0]);
	result.RemoveParticles();

	for(size_t i=0; i<src.NumParticle(); i++)
	{
		GeometryVector OriginalRelative=(src.GetRelativeCoordinates(i))*(1.0/dNumberInEachSide);
		const char * name=src.GetCharacteristics(i).name;

		//d-dimensional loop
		std::vector<size_t> indexes(dim, 0);
		while(indexes.back()!=NumberInEachSide)
		{
			GeometryVector a(OriginalRelative);
			for(DimensionType j=0; j<dim; j++)
				a.x[j]+=indexes[j]/dNumberInEachSide;
			result.Insert(name, a);

			//loop end
			indexes[0]++;
			for(DimensionType j=0; j<dim-1; j++)
			{
				if(indexes[j]==NumberInEachSide)
				{
					indexes[j]=0;
					indexes[j+1]++;
				}
			}
		}
	}
	return result;
}
Configuration MultiplicateStructure(const Configuration & src, const std::vector<size_t> & NumberInEachSide)
{
	assert(NumberInEachSide.size()==src.GetDimension());
	std::vector<double> dNumberInEachSide;
	for(auto iter=NumberInEachSide.begin(); iter!=NumberInEachSide.end(); iter++)
		dNumberInEachSide.push_back((double)(*iter));

	DimensionType dim=src.GetDimension();
	std::vector<GeometryVector> bases;
	for(DimensionType i=0; i<dim; i++)
		bases.push_back(src.GetBasisVector(i)*dNumberInEachSide[i]);
	Configuration result(src);
	result.ChangeBasisVector(&bases[0]);
	result.RemoveParticles();

	for(size_t i=0; i<src.NumParticle(); i++)
	{
		GeometryVector OriginalRelative=(src.GetRelativeCoordinates(i));
		for(int i=0; i< src.GetDimension(); i++)
			OriginalRelative.x[i]/=dNumberInEachSide[i];

		const char * name=src.GetCharacteristics(i).name;

		//d-dimensional loop
		std::vector<size_t> indexes(dim, 0);
		while(indexes.back()!=NumberInEachSide[dim-1])
		{
			GeometryVector a(OriginalRelative);
			for(DimensionType j=0; j<dim; j++)
				a.x[j]+=indexes[j]/dNumberInEachSide[j];
			result.Insert(name, a);

			//loop end
			indexes[0]++;
			for(DimensionType j=0; j<dim-1; j++)
			{
				if(indexes[j]==NumberInEachSide[j])
				{
					indexes[j]=0;
					indexes[j+1]++;
				}
			}
		}
	}
	return result;
}
Configuration GetUnitCubicBox(DimensionType d, double CellSize)
{
	std::vector<GeometryVector> vbas;
	for(DimensionType i=0; i<d; i++)
	{
		GeometryVector t(d);
		t.x[i]=1.0;
		vbas.push_back(t);
	}
	Configuration result(d, &vbas[0], CellSize);

	return result;
}

Configuration GeneratePoisson(DimensionType d, size_t N, RandomGenerator & gen)
{
	Configuration result = GetUnitCubicBox(d);
	result.Resize(N);
	for (size_t i = 0; i < N; i++)
		result.Insert("A", gen);

	return result;
}

#include "DisplaySpheres.h"

double defaultRadius(const Configuration & list)
{
	double typicalLength = std::pow(list.PeriodicVolume() / list.NumParticle(), 1.0 / list.GetDimension());
	double MinLength = 3 * typicalLength;//distance between nearest particles
	for (size_t i = 0; i < list.NumParticle() && i<100; i++)
	{
		list.IterateThroughNeighbors(i, MinLength, [&MinLength, &i](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
		{
			if (SourceAtom != i)
			{
				double l = std::sqrt(shift.Modulus2());
				if (l < MinLength)
					MinLength = l;
			}
		});
	}
	if (MinLength > 0.3*typicalLength && MinLength < typicalLength)
		//we could make sure no spheres overlap (desirable!) by slightly decreasing sphere radius, 
		//so do it!
		typicalLength = MinLength;
	return 0.3*typicalLength;
}

//convert a configuration to a list of spheres and lines
//if BondLength is positive, also draw bond between particles within this distance
//if displacement is not empty, draw an arrow on each sphere to denote this displacement
void ConfigurationToSpheres(const Configuration & list, std::vector<sphere> & spheres, std::vector<line> & lines, double Radius, double BondLength, const std::vector<GeometryVector> & displacement)
{
	GeometryVector DisplayShift;
	double one_over_scale= std::pow(list.PeriodicVolume(), -1.0/list.GetDimension());
	double defaultR = defaultRadius(list);

	if (list.GetDimension() == 3)
	{
		std::vector<GeometryVector> a(8, GeometryVector(0.0, 0.0, 0.0));
		a[1] = a[0] + list.GetBasisVector(0);
		a[2] = a[0] + list.GetBasisVector(1);
		a[3] = a[0] + list.GetBasisVector(2);
		a[4] = a[1] + list.GetBasisVector(1);
		a[5] = a[1] + list.GetBasisVector(2);
		a[6] = a[2] + list.GetBasisVector(2);
		a[7] = a[4] + list.GetBasisVector(2);
		DisplayShift = a[7] * 0.5;
		//one_over_scale = 1.0 / std::sqrt(DisplayShift.Modulus2());
		for (int i = 0; i < 8; i++)
		{
			a[i].MinusFrom(DisplayShift);
			a[i].MultiplyFrom(one_over_scale);
		}
		int point1[12] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6 };
		int point2[12] = { 1, 2, 3, 4, 5, 4, 6, 5, 6, 7, 7, 7 };
		for (int i = 0; i < 12; i++)
		{
			line temp;
			temp.x1 = a[point1[i]].x[0];
			temp.y1 = a[point1[i]].x[1];
			temp.z1 = a[point1[i]].x[2];
			temp.x2 = a[point2[i]].x[0];
			temp.y2 = a[point2[i]].x[1];
			temp.z2 = a[point2[i]].x[2];
			temp.blue = 0.0;
			temp.red = 0.0;
			temp.green = 0.0;
			temp.transparency = 1.0;
			temp.width = 0.5*defaultR*one_over_scale;
			lines.push_back(temp);
		}
	}
	else if (list.GetDimension() == 2)
	{
		std::vector<GeometryVector> a(4, GeometryVector(0.0, 0.0));
		a[1] = a[0] + list.GetBasisVector(0);
		a[2] = a[0] + list.GetBasisVector(1);
		a[3] = a[1] + list.GetBasisVector(1);
		DisplayShift = a[3] * 0.5;
		//one_over_scale = 1.0 / std::sqrt(DisplayShift.Modulus2());
		for (int i = 0; i < 4; i++)
		{
			a[i].MinusFrom(DisplayShift);
			a[i].MultiplyFrom(one_over_scale);
		}
		int point1[4] = { 0, 0, 1, 2 };
		int point2[4] = { 1, 2, 3, 3 };
		for (int i = 0; i < 4; i++)
		{
			line temp;
			temp.x1 = a[point1[i]].x[0];
			temp.y1 = a[point1[i]].x[1];
			temp.z1 = 0.0;
			temp.x2 = a[point2[i]].x[0];
			temp.y2 = a[point2[i]].x[1];
			temp.z2 = 0.0;
			temp.blue = 0.0;
			temp.red = 0.0;
			temp.green = 0.0;
			temp.transparency = 1.0;
			temp.width = 0.5*defaultR*one_over_scale;
			lines.push_back(temp);
		}
	}
	else
	{
		std::cerr << "Error in ConfigurationtoSpheres : unsupported dimension!\n";
		return;
	}

	for (size_t i = 0; i < list.NumParticle(); i++)
	{
		GeometryVector loc = (list.GetCartesianCoordinates(i) - DisplayShift)*one_over_scale;
		sphere temp;
		temp.x = loc.x[0];
		temp.y = loc.x[1];
		temp.z = loc.x[2];

		char chr1 = list.GetCharacteristics(i).name[0];
		if (chr1 >= 'A' && chr1 <= 'Z')
		{
			temp.red = ((chr1 - 'A') / 6)*0.1 + 0.3;
			temp.green = ((chr1 - 'A') % 6)*0.1 + 0.3;
			char chr2 = list.GetCharacteristics(i).name[1];
			if (chr2 >= 'a' && chr2 <= 'z')
				temp.blue = (chr2 - 'a')*0.02 + 0.4;
			else
				temp.blue = 0.35;
			//debug temp
			//if(chr1=='A')
			//{
			//	temp.blue=1.0;
			//	temp.red=1.0;
			//	temp.green=1.0;
			//}
		}
		else
		{
			temp.red = 0.25;
			temp.green = 0.25;
			temp.blue = 0.25;
		}
		temp.transparency = 1.00;

		if (Radius == 0.0)
			temp.radius = list.GetCharacteristics(i).radius*one_over_scale;
		else
			temp.radius = Radius*one_over_scale;

		spheres.push_back(temp);

		if (BondLength > 0)
		{
			auto IterateFunc = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
			{
				//todo: does it plot two lines per bond? If so, correct it.
				if (shift.Modulus2() == 0.0)
					return;
				GeometryVector end = loc + shift*one_over_scale;
				line temp;
				temp.x1 = loc.x[0];
				temp.y1 = loc.x[1];
				temp.z1 = loc.x[2];
				temp.x2 = end.x[0];
				temp.y2 = end.x[1];
				temp.z2 = end.x[2];
				temp.blue = 0.25;
				temp.red = 0.25;
				temp.green = 0.25;
				temp.transparency = 1.0;
				temp.width = 0.2*spheres[0].radius;
				lines.push_back(temp);
			};
			list.IterateThroughNeighbors(i, BondLength, IterateFunc);
		}

		if (i < displacement.size())
		{
			if (displacement[i].Modulus2()>1e-5/one_over_scale/one_over_scale)
			{
				GeometryVector loc2 = loc + displacement[i] * one_over_scale;
				line temp;
				temp.x1 = loc.x[0];
				temp.y1 = loc.x[1];
				temp.z1 = loc.x[2];
				temp.x2 = loc2.x[0];
				temp.y2 = loc2.x[1];
				temp.z2 = loc2.x[2];

				if (list.GetDimension() == 2)
					temp.z1 = temp.z2 = -1.0;

				temp.blue = 1.00;
				temp.red = 1.00;
				temp.green = 0.00;
				temp.transparency = 1.0;
				temp.width = 0.4*spheres[0].radius;
				temp.arrowHead = true;
				lines.push_back(temp);
			}
		}
	}
}

void DisplayConfiguration(const Configuration & a, double BondLength, double Radius)
{
	std::vector<sphere> spheres;
	std::vector<line> lines;

	ConfigurationToSpheres(a, spheres, lines, Radius, BondLength);


	DisplaySpheres(spheres, lines);
}

void Display3DConfigurationMovie(ConfigurationPack & pk, double BondLength, double Radius)
{
	if (pk.NumConfig() == 0)
	{
		std::cout << "warning in Display3DConfigurationMovie : ConfigurationPack is empty. Returning!\n";
		return;
	}
	Configuration a = pk.GetConfig(0);

	auto GetFrameFunc = [&a, &pk, &Radius, &BondLength](std::vector<sphere> & spheres, std::vector<line> & lines, long FrameNumber) ->void
	{
		spheres.clear();
		lines.clear();
		Configuration l = pk.GetConfig(FrameNumber);
		ConfigurationToSpheres(l, spheres, lines, Radius, BondLength);
	};
	
	DisplaySpheres_Movie(GetFrameFunc, pk.NumConfig());
}

#ifdef USE_OPENGL
#include "StructureFactor.h"
#include <GL/freeglut.h>
void Display3DStructureFactor(const Configuration& a, double KCutOff, double KSphereRadius)
{
	assert(a.GetDimension()==3);
	double Smin=1000000, Smax=-1;
	auto SetColorFunc = [&] (double Sk, sphere & s) ->void
	{
		s.red=(Sk-Smin)/(Smax-Smin);
		s.green=0.0;
		s.blue=1.0-s.red;
	};
	auto DrawColorBarFunc = [&] (size_t w, size_t h) ->void
	{
		glBegin(GL_POLYGON);
		glColor3f(0.0, 0.0, 1.0);
		glVertex2f(0.87*w, 0.1*h);
		glVertex2f(0.9*w, 0.1*h);
		glColor3f(1.0, 0.0, 0.0);
		glVertex2f(0.9*w, 0.9*h);
		glVertex2f(0.87*w, 0.9*h);
		glEnd();
		std::stringstream smin;
		smin<<Smin;
		DisplayText(smin.str(), 0.91*w, 0.1*h);
		std::stringstream smax;
		smax<<Smax;
		DisplayText(smax.str(), 0.91*w, 0.9*h);
	};
	GeometryVector recip[3];
	for(int i=0; i<3; i++)
		recip[i]=a.GetReciprocalBasisVector(i);
	PeriodicCellList<Empty> lReciprocal(3, recip, MaxDistance, false);
	lReciprocal.Insert(Empty(), GeometryVector(3));
	std::vector<sphere> VSpheres;
	auto IterateFunc = [&] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle)
	{
		if(shift.Modulus2()==0.0)
			return;
		sphere s;
		s.x=shift.x[0]/KCutOff;
		s.y=shift.x[1]/KCutOff;
		s.z=shift.x[2]/KCutOff;
		double sk=std::log10( StructureFactor(a, shift) );
		if(sk<Smin)
			Smin=sk;
		if(sk>Smax)
			Smax=sk;
		s.transparency=sk;//temperarily use this data member to store S(k)
		VSpheres.push_back(s);
	};
	lReciprocal.IterateThroughNeighbors(GeometryVector(3), KCutOff, IterateFunc);
	if(VSpheres.empty())
		return;
	Smax=std::ceil(Smax);
	Smin=std::floor(Smin);
	if(KSphereRadius==0.0)
	{
		//determine the correct radius of spheres
		auto CompareFunc = [] (const sphere & left, const sphere & right) ->bool
		{
			return left.x*left.x+left.y*left.y+left.z*left.z<right.x*right.x+right.y*right.y+right.z*right.z;
		};
		std::nth_element(VSpheres.begin(), VSpheres.begin(), VSpheres.end(), CompareFunc);
		sphere s=VSpheres[0];
		KSphereRadius = 0.2 * std::sqrt( s.x*s.x+s.y*s.y+s.z*s.z );
	}
	for(auto iter=VSpheres.begin(); iter!=VSpheres.end(); iter++)
	{
		SetColorFunc(iter->transparency, *iter);

		//set transparency accoring to intensity
		//this enables us to view non-zero values better
		iter->transparency=(iter->transparency-Smin)/(Smax-Smin)*0.8+0.2;
		//iter->transparency=0.5;

		iter->radius=KSphereRadius;
	}
	{
		//use a black sphere to indicate origin
		sphere s=VSpheres[0];
		s.x=s.y=s.z=s.red=s.green=s.blue=0.0;
		s.transparency=1.0;
		VSpheres.push_back(s);
	}
	DisplaySpheres(VSpheres, std::vector<line>(), DrawColorBarFunc);
}
#else

void Display3DStructureFactor(const Configuration& a, double KCutOff, double KSphereRadius)
{
	std::cerr<<"Error in Display3DStructureFactor : OpenGL not enabled!\n";
}
#endif

void ConfigurationPack::Open(const std::string & prefix)
{
	this->IndexName=prefix+".ConfigPackIndex";
	this->PackName=prefix+".ConfigPack";
	std::fstream ifile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	if (ifile.good())
	{
		ifile.read((char*)(&NConfig), sizeof(NConfig));
		if (ifile.eof())
			NConfig = 0;
	}
	else
		NConfig=0;
}
long long ConfigurationPack::NumConfig(void) const
{
	return NConfig;
}
void ConfigurationPack::AddConfig(const Configuration & c)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary);
	if(indexFile.good()==false)
	{
		indexFile.clear();
		indexFile.close();
		indexFile.clear();
		indexFile.open(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	}
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::ate);
	if(packFile.good()==false)
	{
		packFile.clear();
		packFile.close();
		packFile.clear();
		packFile.open(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	}
	this->NConfig++;
	indexFile.seekp(0, std::fstream::beg);
	indexFile.write( (char*)(&NConfig), sizeof(NConfig) );
	long long loc=packFile.tellp();
	indexFile.seekp(0, std::fstream::end);
	indexFile.write( (char*)(&loc), sizeof(loc) );
	c.WriteBinary(packFile);
	indexFile.close();
	if (indexFile.fail())
		std::cerr << "Warning in ConfigurationPack::AddConfig : Index file failing!\n";
	packFile.close();
	if (packFile.fail())
		std::cerr << "Warning in ConfigurationPack::AddConfig : Pack file failing!\n";
}
Configuration ConfigurationPack::GetConfig(long long i) const
{
	if (i >= this->NConfig)
	{
		std::cerr << "Error in ConfigurationPack::GetConfig : Index is larger than number of configuration!\n";
		exit(1);
	}
	std::fstream indexFile(IndexName.c_str(), std::fstream::in | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::in | std::fstream::binary);
	long long loc;
	indexFile.seekg(sizeof(long long)*(i + 1));
	indexFile.read((char*)(&loc), sizeof(loc));
	packFile.seekg(std::streampos(loc));

	DimensionType dim;
	packFile.read((char *)(&dim), sizeof(dim));
	if (dim < 1 || dim > ::MaxDimension)
		return Configuration();

	packFile.seekg(std::streampos(loc));

	try
	{
		return packFile.good() ? Configuration(packFile) : Configuration();
	}
	catch (std::exception a)
	{
		std::cerr << a.what();
		return Configuration();
	}
}
void ConfigurationPack::Clear(void)
{
	std::fstream indexFile(IndexName.c_str(), std::fstream::out | std::fstream::binary);
	std::fstream packFile(PackName.c_str(), std::fstream::out | std::fstream::binary | std::fstream::ate);
	this->NConfig=0;
}


//polydisperse sphere packing
SpherePacking::SpherePacking(DimensionType Dimension, GeometryVector * BasisVectors, double CellSize, bool UseSortedList) : PeriodicCellList<double>(Dimension, BasisVectors, CellSize, UseSortedList)
{
	//std::cerr<<"SpherePacking CellSize="<<CellSize<<'\n';
	MaxRadius=0;
}

void SpherePacking::Insert(double radius, const GeometryVector & RelativeCoordinate)
{
	UpdateMaxRadius();
	this->ParticleCartesians.push_back(GeometryVector(this->Dimension));
	this->ParticleRelatives.push_back(GeometryVector(this->Dimension));
	this->ParticleCharacteristics.push_back(radius);
	for(DimensionType i=0; i<this->Dimension; i++)
	{
		assert( RelativeCoordinate.x[i]==RelativeCoordinate.x[i] );

		//NEED TO subtract the floor TWICE
		this->ParticleRelatives.back().x[i]=RelativeCoordinate.x[i]-std::floor(RelativeCoordinate.x[i]);
		this->ParticleRelatives.back().x[i]-=std::floor(this->ParticleRelatives.back().x[i]);
	}
	this->UpdateCartesianCoordinates(this->NumParticle()-1);
	this->Cells[this->GetIndex(this->NumParticle()-1)].push_back(this->NumParticle()-1);

	if(this->MaxRadius<radius)
		this->MaxRadius=radius;
}
//relative coord and halfsize 
//Radius : cartesian radius of the upcoming sphere
bool SpherePacking::CheckVoxelOverlap(const GeometryVector & cr, double halfsize, double Radius) const
{
	UpdateMaxRadius();
	bool Overlap=false;
	const GeometryVector * ba=this->BasisVector;
	DimensionType dim=this->Dimension;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius, [&Overlap, &halfsize, &ba, dim, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
	{
		//additional processing of shift
		GeometryVector ss(shift);
		for(DimensionType d=0; d<dim; d++)
			if(ss.Dot(ba[d])>0.0)
				ss.AddFrom( halfsize*ba[d]);
			else
				ss.MinusFrom( halfsize*ba[d]);

		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);
		if(ss.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	return Overlap;
}

//returns 0 if voxel is completely empty
//2 if voxel is fully occupied
//1 otherwise
int SpherePacking::CheckVoxelOverlap2(const GeometryVector & cr, double halfsize, double Radius) const
{		
	UpdateMaxRadius();
	bool Overlap=false;
	bool PartialOverlap=false;
	const GeometryVector * ba=this->BasisVector;
	DimensionType dim=this->Dimension;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	double maxVoxelDiagonal=this->MaxDiagonalLength()*halfsize;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius+maxVoxelDiagonal, [&maxVoxelDiagonal, &PartialOverlap, &Overlap, &halfsize, &ba, dim, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
	{
		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);

		double temp=OverlapRadius+maxVoxelDiagonal;
		if(shift.Modulus2()<temp*temp)
			PartialOverlap=true;
		//additional processing of shift
		GeometryVector ss(shift);
		for(DimensionType d=0; d<dim; d++)
			if(ss.Dot(ba[d])>0.0)
				ss.AddFrom( halfsize*ba[d]);
			else
				ss.MinusFrom( halfsize*ba[d]);

		if(ss.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	if(Overlap)
		return 2;
	else if(PartialOverlap)
		return 1;
	else
		return 0;
}


//relative coord and halfsize 
//Radius : cartesian radius of the upcoming sphere
bool SpherePacking::CheckOverlap(const GeometryVector & cr, double Radius) const
{
	UpdateMaxRadius();
	bool Overlap=false;
	const std::vector<double> * pSphereRadii = & this->ParticleCharacteristics;
	this->IterateThroughNeighbors(cr, Radius+MaxRadius, [&Overlap, &Radius, &pSphereRadii](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
	{
		double OverlapRadius=Radius+pSphereRadii->at(Sourceparticle);
		if(shift.Modulus2()<OverlapRadius*OverlapRadius)
			Overlap=true;
	}, &Overlap);
	return Overlap;
}


//pixelizing the configuration to MeshSide pixes in each side
void DigitizeConfiguration(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> & occu)
{
	DimensionType dim = c.GetDimension();
	size_t MeshTotal = std::pow(MeshSide, dim);
	struct VicinityMesh
	{
		signed long DeltaIndex[::MaxDimension];
	};

	//get the list of mesh points occupied by a sphere
	std::vector<VicinityMesh> mesh;

	GeometryVector newbas[::MaxDimension];
	for (int i = 0; i<dim; i++)
		newbas[i] = c.GetBasisVector(i)*(1.0 / MeshSide);
	Configuration * pt = new Configuration(dim, newbas, ::MaxDistance, false);
	pt->Insert("A", GeometryVector(dim));
	pt->IterateThroughNeighbors(GeometryVector(dim), radius, [&](const GeometryVector & s, const GeometryVector & l, const signed long * p, size_t source) ->void
	{
		VicinityMesh temp;
		for (int i = 0; i<dim; i++)
			temp.DeltaIndex[i] = p[i];
		mesh.push_back(temp);
	});
	delete pt;

	occu.assign(MeshTotal, 0);
	signed long m = MeshSide;
	int num = c.NumParticle();
	for (int i = 0; i<num; i++)
	{
		GeometryVector rel = c.GetRelativeCoordinates(i);
		for (auto iter = mesh.begin(); iter != mesh.end(); ++iter)
		{
			size_t index = 0;
			for (int j = 0; j<dim; j++)
			{
				index *= MeshSide;
				signed long temp = std::floor(rel.x[j] * MeshSide);
				temp += iter->DeltaIndex[j];
				while (temp >= m)
					temp -= m;
				while (temp<0)
					temp += m;
				index += temp;
			}
			occu[index] = 1;
		}
	}
}
//calculate volume fraction if each particle is replaced with a hypersphere of certain radius
//do this by pixelizing the configuration to MeshSide pixes in each side
double Volume(const Configuration & c, double radius, size_t MeshSide, std::vector<VoxelType> * pbuffer)
{
	size_t MeshTotal = std::pow(MeshSide, c.GetDimension());
	size_t InsideCount=0;
	if (pbuffer == nullptr)
	{
		std::vector<VoxelType> occu;
		DigitizeConfiguration(c, radius, MeshSide, occu);
		for (auto iter = occu.begin(); iter != occu.end(); ++iter)
			InsideCount += *iter;
	}
	else
	{
		std::vector<VoxelType> & occu = *pbuffer;
		DigitizeConfiguration(c, radius, MeshSide, occu);
		for (auto iter = occu.begin(); iter != occu.end(); ++iter)
			InsideCount += *iter;
	}
	return c.PeriodicVolume()*InsideCount/MeshTotal;
}



#include <boost/multi_index_container.hpp>
#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <set>
//return the minimum D such that if each particle is replaced with a sphere with diameter D, particles a and b will be connected.
double MinConnectDiameter(const Configuration & c, size_t a, size_t b)
{
	double clength = 0.0;
	double mlength = 2 * std::pow(c.PeriodicVolume() / c.NumParticle(), 1.0 / c.GetDimension());
	std::set<size_t> cluster;
	struct toSearchStruct
	{
		double dist;
		size_t i;
		toSearchStruct(double dist, size_t i) : dist(dist), i(i)
		{}
	};
	typedef boost::multi_index_container<
		toSearchStruct,
		boost::multi_index::indexed_by<
		boost::multi_index::ordered_unique<boost::multi_index::member<toSearchStruct, size_t, &toSearchStruct::i> >,
		boost::multi_index::ordered_non_unique<boost::multi_index::member<toSearchStruct, double, &toSearchStruct::dist> >
		>
	> toSearch_Set;
	toSearch_Set toAdd;
	toSearch_Set::nth_index<0>::type & iview = toAdd.get<0>();
	toSearch_Set::nth_index<1>::type & dview = toAdd.get<1>();

	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
	{
		if (cluster.find(SourceAtom) == cluster.end())
		{
			double dist = std::sqrt(shift.Modulus2());
			auto iter = iview.find(SourceAtom);
			if (iter == iview.end())
				toAdd.insert(toSearchStruct(dist, SourceAtom));
			else if (iter->dist > dist)
			{
				toSearchStruct temp(dist, SourceAtom);
				toAdd.replace(iter, temp);
			}
		}
	};
	cluster.insert(a);
	GeometryVector rel = c.GetRelativeCoordinates(a);
	c.IterateThroughNeighbors(rel, mlength, IterateFunction);
	for (;;)
	{
		if (toAdd.empty())
		{
			mlength *= 2;
			for (auto iter = cluster.begin(); iter != cluster.end(); iter++)
			{
				GeometryVector rel = c.GetRelativeCoordinates(*iter);
				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
			}
			//std::cout<<"Warning in MinConnectDiameter : resizing mlength. toAdd size="<<toAdd.size()<<'\n';
		}
		else
		{
			auto ibegin = dview.begin();
			toSearchStruct temp = *ibegin;
			dview.erase(ibegin);
			clength = std::max(clength, temp.dist);
			size_t np = temp.i;
			if (np == b)
				return clength;
			if (cluster.find(np) == cluster.end())
			{
				cluster.insert(np);
				GeometryVector rel = c.GetRelativeCoordinates(np);
				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
				//std::cerr<<"Explore particle "<<np<<'\n';
			}
		}
	}

	//shouldn't reach this line
	std::cerr << "Error in MinConnectDiameter\n";
	return 0.0;
}
//double oldMinConnectDiameter(const Configuration & c, size_t a, size_t b)
//{
//	double clength=0.0;
//	double mlength=1.5*std::pow(c.PeriodicVolume()/c.NumParticle(), 1.0/c.GetDimension());
//	std::set<size_t> cluster;
//	std::multimap<double, size_t> toAdd;
//	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
//	{
//		if(cluster.find(SourceAtom)==cluster.end())
//		{
//			double dist=std::sqrt(shift.Modulus2());
//			toAdd.insert(std::make_pair(dist, SourceAtom));
//		}
//	};
//	cluster.insert(a);
//	GeometryVector rel=c.GetRelativeCoordinates(a);
//	c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//	for(;;)
//	{
//		if(toAdd.empty())
//		{
//			mlength*=2;
//			for(auto iter=cluster.begin(); iter!=cluster.end(); iter++)
//			{
//				GeometryVector rel=c.GetRelativeCoordinates(*iter);
//				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//			}
//			//std::cout<<"Warning in MinConnectDiameter : resizing mlength. \n";
//		}
//		else
//		{
//			auto ibegin=toAdd.begin();
//			std::pair<double, size_t> temp=*ibegin;
//			toAdd.erase(ibegin);
//			clength=std::max(clength, temp.first);
//			size_t np=temp.second;
//			if(np==b)
//				return clength;
//			if(cluster.find(np)==cluster.end())
//			{
//				cluster.insert(np);
//				GeometryVector rel=c.GetRelativeCoordinates(np);
//				c.IterateThroughNeighbors(rel, mlength, IterateFunction);
//			}
//		}
//	}
//}

//return the minimum D such that if each particle is replaced with a sphere with diameter D, 
//particle n will connect to its periodic image in direction dir.
double PercolationDiameter(const Configuration & c, size_t n, DimensionType dir)
{
	std::vector<size_t> m(c.GetDimension(), 1);
	m[dir] = 2;
	Configuration t = MultiplicateStructure(c, m);
	return MinConnectDiameter(t, 2 * n, 2 * n + 1);
}

//Replace each particle in c with a hypersphere with radius,
//what's the cluster size InitParticle is in?
//FoundPeriodicImages: whether or not found more than 1 periodic images of a particle, 
//which means the cluster size is actually infinite rather than the returned value.
#include <list>
size_t ClusterSize(const Configuration & c, size_t InitParticle, double radius, bool & FoundPeriodicImages)
{
	FoundPeriodicImages = false;
	std::map<size_t, GeometryVector> Found;
	struct toSearchStruct
	{
		size_t i;
		GeometryVector LatticeShift;
		toSearchStruct(size_t i, const GeometryVector & LatticeShift) : i(i), LatticeShift(LatticeShift)
		{}
	};
	std::list<toSearchStruct> toSearch;
	GeometryVector alreadyLatticeShift = GeometryVector(c.GetDimension());
	auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t SourceAtom)->void
	{
		auto iter = Found.find(SourceAtom);
		GeometryVector l = alreadyLatticeShift + LatticeShift;
		if (iter != Found.end())
		{
			if (!FoundPeriodicImages)
			{
				//this particle has been found before. Test if they belong to the same periodic cell
				GeometryVector deltaLatticeShift = l - iter->second;
				if (deltaLatticeShift.Modulus2() > LengthPrecision*LengthPrecision)
					FoundPeriodicImages = true;
			}
		}
		else
		{
			toSearch.push_back(toSearchStruct(SourceAtom, l));
			Found.insert(std::make_pair(SourceAtom, l));
		}
	};
	Found.insert(std::make_pair(InitParticle, GeometryVector(c.GetDimension())));
	c.IterateThroughNeighbors(InitParticle, radius, IterateFunction);
	while (!toSearch.empty())
	{
		toSearchStruct s = toSearch.front();
		toSearch.pop_front();
		alreadyLatticeShift = s.LatticeShift;
		c.IterateThroughNeighbors(s.i, radius, IterateFunction);
	}
	return Found.size();
}


#include "Voronoi.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
std::complex<double> Psi6(const PeriodicCellList<Empty> & Config, size_t j, int NeighborType, signed int m)
{
	typedef std::complex<double> complex;
	std::vector<GeometryVector> NeighborDisplacements;
	DimensionType dim = Config.GetDimension();
	if (NeighborType == 0)
	{
		std::vector<GeometryVector> neighbors;
		double TypicalLength = std::pow(Config.PeriodicVolume() / Config.NumParticle(), 1.0 / Config.GetDimension())*1.5;
		double l = TypicalLength;
		size_t Nneighbor;
		if (dim == 2)
			Nneighbor = 6;
		else if (dim == 3)
			Nneighbor = 12;
		else
			std::cerr << "Error in Psi6 : unsupported dimension!\n";
		while (neighbors.size() < Nneighbor+1)
		{
			neighbors.clear();
			Config.IterateThroughNeighbors(j, l, [&neighbors](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) ->void
			{
				neighbors.push_back(shift);
			});
			l *= 2;
		}
		auto CompareFunc = [](const GeometryVector & left, const GeometryVector & right) ->bool
		{
			return left.Modulus2() < right.Modulus2();
		};
		std::sort(neighbors.begin(), neighbors.end(), CompareFunc);
		for (int i = 1; i < Nneighbor + 1; i++)
			NeighborDisplacements.push_back(neighbors[i]);
	}
	else if (NeighborType == 1)
	{
		//std::vector<GeometryVector> neighbors;
		//DimensionType Dimension = Config.GetDimension();
		//std::vector<GeometryVector> vertices;
		//std::vector<signed long> Linkage;
		//GetVoronoiCell(Config, j, vertices, &Linkage);
		//for (size_t i = 0; i < vertices.size(); i++)
		//{
		//	for (DimensionType k = 0; k < Dimension + 1; k++)
		//	{
		//		if (Linkage[i*(Dimension + 1) + k] == -1)
		//			continue;
		//		if (Linkage[i*(Dimension + 1) + k] >= i)
		//			continue;
		//		GeometryVector cj = Config.GetCartesianCoordinates(j);
		//		GeometryVector vi = vertices[i];
		//		GeometryVector vk = vertices[Linkage[i*(Dimension + 1) + k]];
		//		GeometryVector cvi = vi - cj;
		//		GeometryVector vik = vk - vi;
		//		GeometryVector ortho = cvi - cvi.Dot(vik) / vik.Modulus2()*vik;
		//		neighbors.push_back(2.0*ortho);
		//	}
		//}

		//debug temp
		//std::cout << "Neighbors contain\n";
		//for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter)
		//	std::cout << *iter;

		//new code
		{
			DimensionType Dimension = Config.GetDimension();
			std::vector<GeometryVector> vertices;
			std::vector<signed long> Linkage;
			std::set<size_t> Neighbors;
			GetVoronoiCell(Config, j, vertices, &Linkage, &Neighbors);

			//debug temp
			//std::cout << "new Neighbors contain\n";
			for (auto iter = Neighbors.begin(); iter != Neighbors.end(); ++iter)
			{
				GeometryVector neighbor = Config.GetRelativeCoordinates(*iter) - Config.GetRelativeCoordinates(j);
				Config.RelativeCoordToMinimumImage(neighbor);
				neighbor = Config.RelativeCoord2CartesianCoord(neighbor);
				NeighborDisplacements.push_back(neighbor);
				//debug temp
				//std::cout << neighbor;
			}
		}

		//complex sum = 0;
		//for (int i = 0; i < neighbors.size(); i++)
		//	sum += std::exp(complex(0.0, 6.0)*std::atan2(neighbors[i].x[1], neighbors[i].x[0]));
		//return sum / (double)(neighbors.size());

	}
	else
	{
		std::cerr << "Error in Psi6 : unsupported NeighborType!\n";
		return complex(0, 0);
	}
	complex sum = 0;
	if (dim == 2)
		for (auto iter = NeighborDisplacements.begin(); iter != NeighborDisplacements.end(); ++iter)
			sum += std::exp(complex(0.0, 6.0)*std::atan2(iter->x[1], iter->x[0]));
	else if (dim == 3)
		for (auto iter = NeighborDisplacements.begin(); iter != NeighborDisplacements.end(); ++iter)
			sum += boost::math::spherical_harmonic(6, m, std::acos(iter->x[2] / std::sqrt(iter->Modulus2())), std::atan2(iter->x[1], iter->x[0]));
			//sum += std::exp(complex(0.0, m)*std::atan2(iter->x[1], iter->x[0])) * LegenderSphericalPlm(6, m, iter->x[2] / std::sqrt(iter->Modulus2()));
	else
		std::cerr << "Error in Psi6 : unsupported dimension!\n";
	return sum / (double)(NeighborDisplacements.size());
}


//display g1(r) of a 3D ensemble
//g1(r) is displayed using spheres in a 3D grid. Sphere volume is proportional to g1.
//sphere color correspond to ln(g1+1), when this quantity increases from 0 to its maximum,
//color changes from black to blue, cyan, yellow, and red
void DisplayG1_3D(ConfigurationPack pk, size_t NGridPerSide)
{
	Configuration list = pk.GetConfig(0);
	if (list.GetDimension() != 3)
	{
		std::cerr << "Error in DisplayG1_3D : Configuration 0 is not three-dimensional.\n";
		return;
	}

	//generate g1(r)
	size_t NBins = NGridPerSide*NGridPerSide*NGridPerSide;
	size_t sumNumParticle = 0;
	std::vector<double> g1(NBins, 0.0);
	for (size_t i = 0; i < pk.NumConfig(); i++)
	{
		Configuration c = pk.GetConfig(i);
		if (c.GetDimension() != 3)
		{
			std::cerr << "Error in DisplayG1_3D : Configuration "<<i<<" is not three-dimensional.\n";
			return;
		}
		sumNumParticle += c.NumParticle();
		for (size_t j = 0; j < c.NumParticle(); j++)
		{
			GeometryVector t = c.GetRelativeCoordinates(j);
			int nx = std::floor(t.x[0] * NGridPerSide);
			int ny = std::floor(t.x[1] * NGridPerSide);
			int nz = std::floor(t.x[2] * NGridPerSide);
			g1[nx*NGridPerSide*NGridPerSide + ny*NGridPerSide + nz] += 1.0;
		}
	}
	double max = 0.0;
	for (auto iter = g1.begin(); iter != g1.end(); ++iter)
	{
		(*iter) *= ((double)(NBins) / sumNumParticle);
		max = std::max(max, *iter);
	}
	std::cout << "max=" << max << '\n';

	//plot g1(r)
	std::vector<line> lines;
	std::vector<sphere> spheres;
	std::vector<GeometryVector> a(8, GeometryVector(0.0, 0.0, 0.0));
	a[1] = a[0] + list.GetBasisVector(0);
	a[2] = a[0] + list.GetBasisVector(1);
	a[3] = a[0] + list.GetBasisVector(2);
	a[4] = a[1] + list.GetBasisVector(1);
	a[5] = a[1] + list.GetBasisVector(2);
	a[6] = a[2] + list.GetBasisVector(2);
	a[7] = a[4] + list.GetBasisVector(2);
	GeometryVector DisplayShift = a[7] * 0.5;
	double one_over_scale = 1.0 / std::sqrt(DisplayShift.Modulus2());
	for (int i = 0; i < 8; i++)
	{
		a[i].MinusFrom(DisplayShift);
		a[i].MultiplyFrom(one_over_scale);
	}
	int point1[12] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6 };
	int point2[12] = { 1, 2, 3, 4, 5, 4, 6, 5, 6, 7, 7, 7 };
	for (int i = 0; i < 12; i++)
	{
		line temp;
		temp.x1 = a[point1[i]].x[0];
		temp.y1 = a[point1[i]].x[1];
		temp.z1 = a[point1[i]].x[2];
		temp.x2 = a[point2[i]].x[0];
		temp.y2 = a[point2[i]].x[1];
		temp.z2 = a[point2[i]].x[2];
		temp.blue = 0.0;
		temp.red = 0.0;
		temp.green = 0.0;
		temp.transparency = 1.0;
		temp.width = 0.3*0.5/NGridPerSide;
		lines.push_back(temp);
	}
	for (int i = 0; i < NGridPerSide; i++)
		for (int j = 0; j < NGridPerSide; j++)
			for (int k = 0; k < NGridPerSide; k++)
			{
				GeometryVector t((0.5 + i) / NGridPerSide, (0.5 + j) / NGridPerSide, (0.5 + k) / NGridPerSide);
				t = list.RelativeCoord2CartesianCoord(t);
				t.MinusFrom(DisplayShift);
				t.MultiplyFrom(one_over_scale);
				double & g = g1[i*NGridPerSide*NGridPerSide + j*NGridPerSide + k];
				sphere temp;
				temp.x = t.x[0];
				temp.y = t.x[1];
				temp.z = t.x[2];
				temp.transparency = 1.0;
				temp.radius = std::pow(g / max, 0.3333)*0.5 / NGridPerSide;
				double tt = std::log(g + 1.0) / std::log(max + 1.0);
				if (tt < 0.25)
				{
					temp.red = 0.0;
					temp.green = 0.0;
					temp.blue = tt / 0.25;
				}
				else if (tt < 0.5)
				{
					temp.red = 0.0;
					temp.green = (tt - 0.25) / 0.25;
					temp.blue = 1.0;
				}
				else if (tt < 0.75)
				{
					temp.red = (tt - 0.5) / 0.25;
					temp.green = 1.0;
					temp.blue = 1.0 - (tt - 0.5) / 0.25;
				}
				else
				{
					temp.red = 1.0;
					temp.green = (tt - 0.75) / 0.25;
					temp.blue = 0.0;
				}
				spheres.push_back(temp);
			}

	DisplaySpheres(spheres, lines);
	return;
}

#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/vector.hpp> 
#include <boost/numeric/ublas/io.hpp> 
#include <boost/numeric/ublas/vector_proxy.hpp> 
#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/triangular.hpp> 
#include <boost/numeric/ublas/lu.hpp> 
namespace
{
	struct ParticlePair
	{
		size_t i, j;
		GeometryVector rij;
		ParticlePair(size_t i, size_t j, const GeometryVector & rij) : i(i), j(j), rij(rij)
		{}
	};

	//make sure that in a particle pair, i's nearest neighbor (other than j) is closer than j's nearest neighbor (other than i)
	void StandardizeParticlePairDirection(const Configuration & c, ParticlePair & src)
	{
		size_t i = src.i, j = src.j;

		for (double factor = 1.0;; factor *= 1.5)
		{
			double maxDist = std::pow(c.PeriodicVolume() / c.NumParticle(), 1.0 / c.GetDimension())*factor;

			//generate the neighbor list of two particles
			std::vector<double> l2s, l2s2;
			auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
			{
				if (Sourceparticle != i && Sourceparticle != j)
					l2s.push_back(shift.Modulus2());
			};
			c.IterateThroughNeighbors(j, maxDist, IterateFunction);
			std::sort(l2s.begin(), l2s.end());
			std::swap(l2s2, l2s);
			c.IterateThroughNeighbors(i, maxDist, IterateFunction);
			std::sort(l2s.begin(), l2s.end());

			//compare the two neighbor lists
			for (size_t i = 0;; i++)
			{
				//when two neighbor lists reach the end simultaneously, we cannot decide an order of the particles
				if (i == l2s.size() && i == l2s2.size())
					break;
				//when one neighbor list reach the end, the other particle comes first
				if (i == l2s.size())
				{
					std::swap(src.i, src.j);
					src.rij.MultiplyFrom(-1.0);
					return;
				}
				else if (i == l2s2.size())
					return;

				if (l2s[i] < l2s2[i])
					return;
				else if (l2s[i]>l2s2[i])
				{
					std::swap(src.i, src.j);
					src.rij.MultiplyFrom(-1.0);
					return;
				}
			}
			if (factor > 5.0)
			{
				std::cerr << "Error in StandardizeParticlePairDirection : unable to determine particle order!\n";
				return;
			}
		}
	}

	//find num shortest particle pairs
	std::vector<ParticlePair> findNearestParticlePairs(const Configuration & c, size_t num)
	{
		if (c.NumParticle() < 2 || c.GetDimension() == 0)
		{
			std::cerr << "Error in findNearestParticlePairs : dimension or number of particle is too small.\n";
			return std::vector<ParticlePair>();
		}

		std::vector<ParticlePair> result;

		for (double factor = 1.0; result.size() < num; factor *= 1.5)
		{
			std::multimap<double, ParticlePair> pairs;
			double maxDist = std::pow(c.PeriodicVolume() / c.NumParticle(), 1.0 / c.GetDimension())*factor;
			for (size_t i = 1; i < c.NumParticle(); i++)
			{
				struct neighbor
				{
					size_t j;
					GeometryVector shift;
					double l2;
					neighbor(size_t j, double l2, GeometryVector shift) : j(j), l2(l2), shift(shift)
					{}
				};
				std::vector<neighbor> neighbors;
				auto IterateFunction = [&](const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t Sourceparticle) -> void
				{
					double l2 = shift.Modulus2();
					if (Sourceparticle<i)
						neighbors.push_back(neighbor(Sourceparticle, l2, shift));
				};
				c.IterateThroughNeighbors(i, maxDist, IterateFunction);
				for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter)
					pairs.insert(std::make_pair(iter->l2, ParticlePair(i, iter->j, iter->shift)));

				while (pairs.size()>num)
				{
					auto iter = pairs.end();
					iter--;
					pairs.erase(iter);
				}

				if (pairs.size() == num)
					maxDist = std::sqrt(pairs.rbegin()->first);
			}

			for (auto iter = pairs.begin(); iter != pairs.end(); ++iter)
				result.push_back(iter->second);
		}

		return result;
	}

	template<class T>
	bool InvertMatrix(const boost::numeric::ublas::matrix<T>& input, boost::numeric::ublas::matrix<T>& inverse)
	{
		using namespace boost::numeric::ublas;
		typedef permutation_matrix<std::size_t> pmatrix;
		// create a working copy of the input 
		matrix<T> A(input);
		// create a permutation matrix for the LU-factorization 
		pmatrix pm(A.size1());
		// perform LU-factorization 
		int res = lu_factorize(A, pm);
		if (res != 0)
			return false;
		// create identity matrix of "inverse" 
		inverse.assign(boost::numeric::ublas::identity_matrix<T>(A.size1()));
		// backsubstitute to get the inverse 
		lu_substitute(A, pm, inverse);
		return true;
	}

	//output a matrix from a configuration
	//if a configuration is rotated or reflected, it is guaranteed that the column vectors of the outputed matrix
	//is simply the rotated or reflected column vectors of the outputed matrix of the original configuration
	//pairs simply outputs the pairs used to construct this matrix
	boost::numeric::ublas::matrix<double> StandardizedNearestNeighborMatrix(const Configuration & c, std::vector<ParticlePair> & pairs)
	{
		boost::numeric::ublas::matrix<double> result(c.GetDimension(), c.GetDimension());
		pairs = findNearestParticlePairs(c, c.GetDimension());
		for (int i = 0; i < c.GetDimension(); i++)
		{
			StandardizeParticlePairDirection(c, pairs[i]);
			GeometryVector r = pairs[i].rij;
			for (int j = 0; j < c.GetDimension(); j++)
				result(j, i) = r.x[j];
		}
		return result;
	}

	//test if there is a translation vector and a rotation/inversion matrix, that every particle in a, after rotation and translation, can find a correspondence in b.
	bool CanFindCorresponding(const Configuration & a, const Configuration & b)
	{
		if (a.GetDimension() != b.GetDimension())
			return false;
		DimensionType d = a.GetDimension();
		size_t N = a.NumParticle();

		//find rotation/inversion matrix
		std::vector<ParticlePair> vp1, vp2;
		boost::numeric::ublas::matrix<double> mr1 = StandardizedNearestNeighborMatrix(a, vp1), mr2 = StandardizedNearestNeighborMatrix(b, vp2);
		boost::numeric::ublas::matrix<double> mT(d, d), temp(d, d);
		mT.clear();
		temp.clear();

		InvertMatrix(mr1, temp);
		boost::numeric::ublas::axpy_prod(mr2, temp, mT, true);

		//debug temp
		//std::cout << mT(0, 0) << " " << mT(0, 1) << " " << mT(1, 0) << " " << mT(1, 1) << '\n';

		//find translation
		auto GetBoostVectorFunc = [&](const GeometryVector & src) -> boost::numeric::ublas::vector<double>
		{
			boost::numeric::ublas::vector<double> result(d);
			for (DimensionType i = 0; i < d; i++)
				result(i) = src.x[i];
			return result;
		};
		boost::numeric::ublas::vector<double> v1 = GetBoostVectorFunc(a.GetCartesianCoordinates(vp1[0].i));
		boost::numeric::ublas::vector<double> v2 = GetBoostVectorFunc(b.GetCartesianCoordinates(vp2[0].i));
		boost::numeric::ublas::vector<double> vtemp(d), Translation(d);
		boost::numeric::ublas::axpy_prod(mT, v1, vtemp, true);
		Translation = v2 - vtemp;

		for (size_t i = 0; i < N; i++)
		{
			//find the coordinate corresponding to particle i in configuration a
			boost::numeric::ublas::vector<double> ci1 = GetBoostVectorFunc(a.GetCartesianCoordinates(i));
			boost::numeric::ublas::axpy_prod(mT, ci1, vtemp, true);
			vtemp = vtemp + Translation;
			//convert to GeometryVector
			GeometryVector ci2(d);
			for (DimensionType j = 0; j < d; j++)
				ci2.x[j] = vtemp(j);
			//find nearest neighbor in configuration b
			GeometryVector ri2 = b.CartesianCoord2RelativeCoord(ci2);
			double dist = b.NearestParticleDistance(ri2);
			if (dist > 1e-4)
			{
				std::cout << dist << '\n';
				return false;
			}
		}

		return true;
	}
};

//test if two configurations are equivalent
//works by identifying the closest pair, second closest pair, ..., dth closest pair,
//and use these pairs to find rotation/inversion matrix and the translation vector.
//Therefore, this function works well only if the configurations are disordered.
bool ConfigurationEquivalent(const Configuration & a, const Configuration & b)
{
	return CanFindCorresponding(a, b) && CanFindCorresponding(b, a);
}
