//#ifndef __Boltzmann_Weight__prob__
//#define __Boltzmann_Weight__prob__

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

//#endif


class Probability
{
	std::vector<double> prob;
public : 
	std::vector<double> Prob();
	void ReadProb( std::string & filename );
};

class Bins
{
	std::vector<int> bin;
public:
	std::vector<int> Bin();
	void ReadBin( std::string & filename );
};


class Dats
{
	std::vector< std::vector<double> > dat;
public:
	std::vector< std::vector<double> > Dat();	
	void ReadDat( std::string & filename );
};
