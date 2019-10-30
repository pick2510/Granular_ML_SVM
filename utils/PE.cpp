#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <Eigen/Dense>
#include "omp.h"
#include "utils.h"


double get_pe_sum(const std::string &filename, const std::map<int, double> &radius){
	std::string line;
	std::ifstream contactfile(filename);
	int num_contacts=0, num_properties = 0;
	if (contactfile.is_open()){
		while (std::getline(contactfile, line)){
			if (num_contacts == 0){
				std::vector<std::string> splitted_line;
				split(line, splitted_line);
				num_properties = splitted_line.size();

			}
			num_contacts++;	
		}
		contactfile.clear();
		contactfile.seekg(0, std::ios::beg);
		Eigen::MatrixXd contact_m(num_contacts, num_properties);



	}
	else {
		std::cout << "ERROR Contactfile " << filename << " is not there";
		std::exit(-1);
	}

}




std::map<int, double> get_radius(const std::string &filename)
{
	std::string line;
	std::ifstream radiusf(filename);
	std::map<int, double> radmap;
	if (radiusf.is_open())
	{
		while (std::getline(radiusf, line))
		{
			int index;
			double val;
			std::vector<std::string> splitted_line;
			split(line, splitted_line);
			std::stringstream(splitted_line[0]) >> index;
			std::stringstream(splitted_line[1]) >> val;
			radmap[index] = val;
		}
	} else {
		std::cout << "NO RADIUSFILE\n";
		std::exit(-1); 
	}
	return radmap;
}



int main(int argc, char *argv[])
{
	std::string chainprefix = "/mnt/DEMDAA/AR=2.56/ForBetweenness4-long/DEM/postchain/", radius_prefix = "/home/strebdom/forOmid/";
	int start, increment, stop;
	std::map<int, double> PE_sum;
	std::cin >> start >> increment >> stop;
	std::stringstream radiusfilename;
	radiusfilename << radius_prefix << "RadiusSorted.txt";
	auto radius = get_radius(radiusfilename.str());
	for (int i = start; i < stop; i += increment)
	{
		std::stringstream fname;
		fname << chainprefix << "contact" << i << ".txt";
	}
	return 0;
}