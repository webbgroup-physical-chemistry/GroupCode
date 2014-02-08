#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <stdio.h>
#include <iomanip>
#include <thread>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif


typedef std::vector< std::vector<double> > Matrix;
typedef std::vector< std::vector<std::string> > DatList;

struct experiment
{
    Matrix dat;
    std::vector<double> prob;
    std::vector<int> bin;
};

struct average
{
    std::vector<double> avg;
    std::vector<double> std;
};

class Probability
{
	std::vector<double> prob;
public : 
	std::vector<double> Prob();
	void ReadProb( std::string & filename );
    void clear();
};

class Bins
{
	std::vector<int> bin;
public:
	std::vector<int> Bin();
	void ReadBin( std::string & filename );
    void clear();
};


class Dats
{
	Matrix dat;
public:
	Matrix Dat();	
	void ReadDat( std::string & filename );
    void clear();
};

class Lists
{
    DatList list;
    experiment frame;
public:
    experiment ListDat();
    void ReadList( std::string & filename );
    void ReadListInputs();
    void FrameProb(Probability probability);
};

class Stats
{
    average boltzmannAverage;
//    std::vector<average> results;
    average * results;
    bool angular;
    experiment frames;
    average bootstrapResults;
public:
    void Average();
    average angleAverage(experiment frames);
    average linearAverage(experiment frames);
    average BoltzmannAverage();
    average BootstrapAverage();
    void resampleAverage(int slot);
    void bootstrap(int nresamples = 10000);
    void bootstrap_threading(int nresamples = 10000);
    void addFrame(experiment addframe,bool angle = false);
    void loop(int stepsize, int stop, int thread);
};  

class WriteOutputs
{
    std::string outfile;
public : 
    void extension( std::string basename, std::string suffix, std::string & outname);
    void backup( std::string filename );
    bool fexists( std::string filename );
    void write(std::string outfile, std::string probfile, average boltzmann, average bootstrap,int nsteps);

};
void foo(int i, int th);
