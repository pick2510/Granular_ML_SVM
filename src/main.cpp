#include "Interfaces.h"

#include <omp.h>
#include <ctime>
#include <iostream>

#include <nlopt.h>

int cli()
{
	char tempstring[1000];
	std::istream & ifile = std::cin;
	std::ostream & ofile = std::cout;
	double BondLength = 0.0, Radius = 0.0;
	for (;;)
	{
		ifile >> tempstring;
		if (strcmp(tempstring, "Exit") == 0)
		{
			return 0;
		}
		else if (strcmp(tempstring, "Verbosity") == 0)
		{
			ifile >> Verbosity;
		}
		else if (strcmp(tempstring, "AlsoWriteEPS") == 0)
		{
			ifile >> AlsoWriteEPS;
		}
		else if (strcmp(tempstring, "TimeLimit") == 0)
		{
			size_t limit;
			ifile >> limit;
			::TimeLimit = ::ProgramStart + limit;
		}
		else if (strcmp(tempstring, "InvStatMech") == 0)
		{
			return InvStatMechCLI();
		}
		else if (strcmp(tempstring, "PairStatistics") == 0)
		{
			return PairStatisticsCLI();
		}
		else if (strcmp(tempstring, "CollectiveCoordinate") == 0)
		{
			return CollectiveCoordinateCLI();
		}
		else if (strcmp(tempstring, "Debug") == 0)
		{
			Debug();
		}
		else if (strcmp(tempstring, "BondLength") == 0)
		{
			ifile >> BondLength;
		}
		else if (strcmp(tempstring, "Radius") == 0)
		{
			ifile >> Radius;
		}
		else if (strcmp(tempstring, "Draw") == 0)
		{
			std::string prefix;
			std::cin >> prefix;
			Configuration c = ReadPos(prefix);

			//if .pos file not found, try to find .ConfigPack file
			if (c.GetDimension() == 0)
			{
				ConfigurationPack pk(prefix);
				if (pk.NumConfig() > 0)
					c = pk.GetConfig(0);
			}
			c.StandardizeOrientation();
			c.Resize(1.0);
			::AlsoWriteEPS = true;
			::Plot(prefix, c);
		}
		else if (strcmp(tempstring, "Display") == 0)
		{
			std::string prefix;
			std::cin >> prefix;
			Configuration c = ReadPos(prefix);
			DisplayConfiguration(c, BondLength, Radius);
		}
		else if (strcmp(tempstring, "DisplayMultiple") == 0)
		{
			std::string prefix;
			std::cin >> prefix;
			Configuration c = ReadPos(prefix);
			DisplayConfiguration(MultiplicateStructure(c, 5), BondLength, Radius);
		}
		else if (strcmp(tempstring, "DisplaySk") == 0)
		{
			std::string prefix;
			std::cin >> prefix;
			Configuration c = ReadPos(prefix);
			double MaxK;
			std::cin >> MaxK;
			Display3DStructureFactor(c, MaxK);
		}
		else if (strcmp(tempstring, "Movie") == 0)
		{
			std::string prefix;
			std::cin >> prefix;
			ConfigurationPack pk(prefix);
			Display3DConfigurationMovie(pk, BondLength, Radius);
		}
		else if (strcmp(tempstring, "Separate") == 0)
		{
			std::string prefix1, prefix2;
			std::cout << "Input ConfigPack=";
			std::cin >> prefix1;
			std::cout << "Output ConfigPack=";
			std::cin >> prefix2;
			size_t m, n;
			std::cout << "number of sub-config per configuration=";
			std::cin >> m;
			std::cout << "the sub-config to extract=";
			std::cin >> n;
			ConfigurationPack pk1(prefix1);
			ConfigurationPack pk2(prefix2);
			pk2.Clear();
			for (int i = 0; i < pk1.NumConfig(); i++)
			{
				Configuration c = pk1.GetConfig(i);
				c = SeparateList(c, m, n);
				pk2.AddConfig(c);
			}
		}
		else if (strcmp(tempstring, "CopyConfigs") == 0)
		{
			std::string prefix1, prefix2;
			std::cout << "Input ConfigPack=";
			std::cin >> prefix1;
			std::cout << "Output ConfigPack=";
			std::cin >> prefix2;
			ConfigurationPack pk1(prefix1);
			ConfigurationPack pk2(prefix2);
			for (int i = 0; i < pk1.NumConfig(); i++)
				pk2.AddConfig(pk1.GetConfig(i));
		}
		else if (strcmp(tempstring, "ConfigPackToNameAndCoordinates") == 0)
		{
			std::string prefix1, prefix2;
			std::cout << "Input ConfigPack=";
			std::cin >> prefix1;
			std::cout << "Output prefix=";
			std::cin >> prefix2;
			ConfigurationPack pk1(prefix1);
			progress_display pd(pk1.NumConfig());
			for (int i = 0; i < pk1.NumConfig(); i++)
			{
				std::stringstream ss;
				ss << prefix2 << i << ".txt";
				std::fstream ofile(ss.str(), std::fstream::out);
				ofile.precision(17);
				Configuration c = pk1.GetConfig(i);
				for (int j = 0; j < c.NumParticle(); j++)
				{
					ofile << c.GetCharacteristics(j);
					for(int k=0; k<c.GetDimension(); k++)
						ofile<< " " << c.GetCartesianCoordinates(j).x[k];
					ofile << std::endl;
				}
				pd++;
			}
		}
		else if (strcmp(tempstring, "RefineBasis") == 0)
		{
			std::string prefix1, prefix2;
			std::cout << "Input ConfigPack=";
			std::cin >> prefix1;
			std::cout << "Output ConfigPack=";
			std::cin >> prefix2;
			ConfigurationPack pk1(prefix1);
			ConfigurationPack pk2(prefix2);
			for (int i = 0; i < pk1.NumConfig(); i++)
			{
				Configuration c = pk1.GetConfig(i);
				c.TryRefineBasisVectors();
				pk2.AddConfig(c);
			}
		}
		else if (strcmp(tempstring, "AddConfigurationToPack") == 0)
		{
			std::string prefix1, prefix2;
			std::cout << "Input Configuration=";
			std::cin >> prefix1;
			std::cout << "Output ConfigPack=";
			std::cin >> prefix2;
			std::ifstream ifile(prefix1 + ".configuration", std::fstream::binary);
			Configuration c(ifile);
			ConfigurationPack pk2(prefix2);
			pk2.AddConfig(c);
		}
		else if (strcmp(tempstring, "SeparateAll") == 0)
		{
			std::string prefix1, prefix2;
			std::cout << "Input ConfigPack=";
			std::cin >> prefix1;
			std::cout << "Output ConfigPack=";
			std::cin >> prefix2;
			size_t m;
			std::cout << "number of sub-config per configuration=";
			std::cin >> m;
			ConfigurationPack pk1(prefix1);
			ConfigurationPack pk2(prefix2);
			pk2.Clear();
			for (int i = 0; i < pk1.NumConfig(); i++)
			{
				Configuration c = pk1.GetConfig(i);
				for (int j = 0; j < m; j++)
				{
					Configuration cc = SeparateList(c, m, j);
					pk2.AddConfig(cc);
				}
			}
		}
		else if (strcmp(tempstring, "PlotG1") == 0)
		{
			std::cout << "Input ConfigPack name, and number of grids per side:";
			std::string name;
			std::cin >> name;
			ConfigurationPack pk(name);
			std::vector<Configuration> a;
			for (int i = 0; i < pk.NumConfig(); i++)
				a.push_back(pk.GetConfig(i));
			size_t n;
			std::cin >> n;
			PlotG1_2D(a, n, name + "_g1", "");
			return 0;
		}
		else if (strcmp(tempstring, "DisplayG1") == 0)
		{
			std::cout << "Input ConfigPack name, and number of grids per side:";
			std::string name;
			std::cin >> name;
			ConfigurationPack pk(name);
			size_t n;
			std::cin >> n;
			DisplayG1_3D(pk, n);
			return 0;
		}
		else if (strcmp(tempstring, "CalculatePrincipalRelaxationTime") == 0)
		{
			std::string name;
			double pmax, pmin;
			std::cin >> name >> pmax >> pmin;
			std::vector< std::vector<GeometryVector> > data;
			ReadGraceData(data, name);
			std::vector< std::vector<double> > x(2, std::vector<double>());
			std::vector<double> y;
			for (auto iter = data[0].begin(); iter != data[0].end(); ++iter)
			{
				if (iter->x[1] > pmin && iter->x[1] < pmax)
				{
					x[0].push_back(1.0);
					x[1].push_back(iter->x[0]);
					y.push_back(std::log(iter->x[1]));
				}
			}
			std::vector<double> c = MultiVariableLinearFit(y, x);
			std::cout << (-1.0) / c[1] << " \t";
		}

		else
			std::cerr << "Unrecognized command!\n";
	}
}


int main(int argc, char ** argv)
{
	//return testSpeed();
	//nlopt_srand(9999999);
	nlopt_srand(99999);

	auto start = std::time(nullptr);
	::ProgramStart=start;
	::TimeLimit=::ProgramStart+86400000;//default time limit of 1000 days
	omp_set_nested(true);
	std::cerr.precision(17);
	std::cout.precision(17);

	struct tm * now = localtime( & start );
	logfile << "============================================================\n"
		<< (now->tm_year + 1900) << '-' 
		<< (now->tm_mon + 1) << '-'
		<<  now->tm_mday << ' '
		<< now->tm_hour <<':'
		<< now->tm_min <<':'
		<< now->tm_sec 
		<< '\n';


	try
	{
		if(argc>1 && strcmp(argv[1], "-debug")==0)
			Debug();
		else
		{
			cli();
		}

	}
	catch(char * error)
	{
		std::cout<<error;
	}
	std::cout<<"time spent: "<<std::time(nullptr)-start<<" seconds\n";
	logfile<<"time spent:"<<std::time(nullptr)-start<< " seconds\n============================================================\n";
	return 0;
	//std::cin.get();
}

