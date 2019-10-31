#include "bits/stdc++.h"
#include <Eigen/Dense>
#include "omp.h"
#include "utils.h"

unsigned int FileRead(std::istream &is, std::vector<char> &buff)
{
	is.read(&buff[0], buff.size());
	return is.gcount();
}

unsigned int CountLines(const std::vector<char> &buff, int sz)
{
	int newlines = 0;
	const char *p = &buff[0];
	for (int i = 0; i < sz; i++)
	{
		if (p[i] == '\n')
		{
			newlines++;
		}
	}
	return newlines;
}

double get_pe_sum(const std::string &filename, const std::map<int, double> &radius)
{
	std::string line;
	std::ifstream contactfile(filename);
	int num_contacts = 0, num_properties = 0;
	const double Ystar = 6.5E11 / (2 * (1 - std::pow(0.25, 2)));
	const double Gstar = 6.5E11 / (4 * (2 - .25) * (1 + .25));

	if (contactfile.is_open())
	{
		std::getline(contactfile, line);
		std::vector<std::string> splitted_line;
		split(line, splitted_line);
		num_properties = splitted_line.size();
		contactfile.clear();
		contactfile.seekg(0, std::ios::beg);
		const int SZ = 1024 * 1024;
		std::vector<char> buff(SZ);
		while (int cc = FileRead(contactfile, buff))
		{
			num_contacts += CountLines(buff, cc);
		}
	}
	else
	{
		std::cout << "ERROR Contactfile " << filename << " is not there";
		std::exit(-1);
	}
	Eigen::MatrixXd contact_m = Eigen::MatrixXd::Constant(num_contacts, num_properties, 0.0);
	contactfile.clear();
	contactfile.seekg(0, std::ios::beg);
	int i = 0;
	while (std::getline(contactfile, line))
	{
		std::stringstream ll(line);
		std::string splitted;
		int j = 0;
		while (std::getline(ll, splitted, ' '))
		{
			std::stringstream(splitted) >> contact_m(i, j);
			j++;
		}
		i++;
	}

	Eigen::MatrixXd calc(contact_m.rows(), 10);
	double Pe_tot{0.0};
	for (int row = 0; row < contact_m.rows(); row++)
	{
		calc(row, calc_layout::rad1) = radius.at(contact_m(row, ContactTXTColumns::p1_id));
		calc(row, calc_layout::rad2) = radius.at(contact_m(row, ContactTXTColumns::p2_id));
		calc(row, calc_layout::Radstar) = (calc(row, calc_layout::rad1) * calc(row, calc_layout::rad2) /
										   (calc(row, calc_layout::rad1) + calc(row, calc_layout::rad2)));
		calc(row, calc_layout::kn) =
			4 / 3.0 * Ystar *
			std::sqrt(calc(row, calc_layout::Radstar) *
					  contact_m(row, ContactTXTColumns::contact_overlap));
		calc(row, calc_layout::kt) =
			8.0 * Gstar *
			std::sqrt(calc(row, calc_layout::Radstar) *
					  contact_m(row, ContactTXTColumns::contact_overlap));
		calc(row, calc_layout::Fnor) =
			std::sqrt(std::pow(contact_m(row, ContactTXTColumns::cn_force_x), 2) +
					  std::pow(contact_m(row, ContactTXTColumns::cn_force_y), 2) +
					  std::pow(contact_m(row, ContactTXTColumns::cn_force_z), 2));
		calc(row, calc_layout::Ftan) =
			std::sqrt(std::pow(contact_m(row, ContactTXTColumns::ct_force_x), 2) +
					  std::pow(contact_m(row, ContactTXTColumns::ct_force_y), 2) +
					  std::pow(contact_m(row, ContactTXTColumns::ct_force_z), 2));
		calc(row, calc_layout::Penor) =
			0.5 * (std::pow(calc(row, calc_layout::Fnor), 2) /
				   calc(row, calc_layout::kn));
		calc(row, calc_layout::Petan) =
			0.5 * (std::pow(calc(row, calc_layout::Ftan), 2) /
				   calc(row, calc_layout::kt));
		calc(row, calc_layout::Petot) = calc(row, calc_layout::Penor) + calc(row, calc_layout::Petan);
		Pe_tot += calc(row, calc_layout::Petot);
	}
	return Pe_tot;
}

void output(std::ofstream &out, std::map<int, double> &PE)
{
	for (auto &elem : PE)
	{
		out << elem.first << " " << elem.second << "\n";
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
	}
	else
	{
		std::cout << "NO RADIUSFILE\n";
		std::exit(-1);
	}
	return radmap;
}

int main(int argc, char *argv[])
{
	std::string chainprefix = "/mnt/DEMDAA/AR=2.56/ForBetweenness4-long/DEM/postchain/", radius_prefix = "/home/strebdom/forOmid/";
	int start, increment, stop;
	std::cin >> start >> increment >> stop;
	std::stringstream radiusfilename;
	std::map<int, double> PE_sum;
	radiusfilename << radius_prefix << "RadiusSorted.txt";
	auto radius = get_radius(radiusfilename.str());
	int round = 0;
#pragma omp parallel for
	for (int i = start; i < stop; i += increment)
	{
		std::stringstream fname;
		fname << chainprefix << "contact" << i << ".txt";
		double ij = get_pe_sum(fname.str(), radius);
#pragma omp critical
		{
			round++;
			PE_sum[i] = ij;
			std::cout << std::setprecision(16);
			std::cout << "i: " << i << " PE: " << ij << " Round: " << round << "\n";
		}
	}
	std::ofstream PEout("PE.txt");
	if (PEout.is_open())
	{
		output(PEout, PE_sum);
	}

	return 0;
}