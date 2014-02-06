#include "prob.h"
#include <stdlib.h>

std::vector<double> Probability::Prob() {return prob;};

void Probability::ReadProb( std::string & filename )
{
	std::cout << "Reading " << filename << "..." << std::endl;
	std::string line;
	double value;
	double sum = 0;
	std::ifstream file (filename.c_str());
	if (file.is_open() )
	{
		while (file.good())
		{
			getline( file, line );
			if (not line.empty())
			{
				std::stringstream linestream( line );
				linestream >> value ;
				prob.push_back(value);
				sum += value;
			}
		}
	}
	else
	{
		std::cout << filename << " does not exist..." << std::endl;
		exit(1);
	}	
	double inverseSum = 1.0/sum;
	for (int i=0; i<prob.size() ; i++)
	{
		prob[i] *= inverseSum;
	}
	return;
};

std::vector<int> Bins::Bin() {return bin;}

void Bins::ReadBin( std::string & filename )
{
	std::cout << "Reading " << filename << "..." << std::endl;
	std::string line;
	int value;
	std::ifstream file (filename.c_str());
	if (file.is_open())
	{
		while (file.good())
		{
			getline( file, line);
			if (not line.empty()) 
			{
				std::stringstream linestream(line);
				linestream >> value;
				bin.push_back(value);
			}
		}
	}
	else
	{
		std::cout << filename << " does not exist..." << std::endl;
		exit(1);
	}
	return;
}

std::vector< std::vector<double> > Dats::Dat(){return dat;}

void Dats::ReadDat( std::string & filename ) 
{
	std::cout << "Reading " << filename << "..." << std::endl;
	std::string line;
	std::vector<double> values;
	std::ifstream file (filename.c_str());
	double n;
	if (file.is_open())
	{
		while (file.good())
		{
			getline( file, line);
			if ( not line.empty())
			{
				std::stringstream linestream(line);
				while (linestream >> n) 
				{
					values.push_back(n);	
				}
				dat.push_back(values);
				values.clear();
			}
		}
	}
	else
	{
		std::cout << filename << " does not exist..." << std::endl;
		exit(1);
	}
	//std::cout << dat[0].size() << std::endl;
	//std::cout << dat.size() << std::endl;
	return;
};

