#include "MyOptionParser.h"
#include "prob.h"

using namespace std; 

int main(int argc, char * argv[])
{
	OptionParser parser;
	parser.assign_header("\
	This is intended to calculate a Boltzmann weigghted average\n\
	given some value over instances, with each instance assigned\n\
	either a bin or a probability.  If each instance is assigned\n\
	a bin, a file governing the probability of each bin is also\n\
	required");
	parser.read_cmdline(argv, argv+argc);
	parser.add_option("-p",\
					  "Probability file, probably a .file.prob",\
                      "test/probability.prob");  
	parser.add_option("-l",\
					  "List of files: <.bin> <.xvg>",\
					  "test/list.file");
    parser.add_booloption("-a",\
                          "Use angular average (periodic)",\
                          false);
    parser.add_option("-o",\
                      "output file name",\
                      "out.boltzmann");
    parser.add_option("-n",\
                      "Number of bootstrapping experiments",\
                      "100");
    parser.add_booloption("-b",\
                          "Do bootstrapping",\
                          false);

                       
	parser.parse_args();
	std::string probfile = parser.Options()[0];
	std::string listfile = parser.Options()[1];
    bool angular = parser.boolOptions()[0];
    std::string outfile = parser.Options()[2];
    std::string donsteps = parser.Options()[3].c_str();
    int nsteps ; 
    std::stringstream linestream(donsteps);
    linestream >> nsteps;
    bool doBootstrapping = parser.boolOptions()[1];

	Probability probability;
	probability.ReadProb( probfile );

    Lists list;
    list.ReadList(listfile);
    list.ReadListInputs();
    list.FrameProb(probability);
    experiment results = list.ListDat();

    Stats stats;
    stats.addFrame(results,angular);
    average bw, bs;
    bw = stats.BoltzmannAverage();
    if (doBootstrapping)
    {
        stats.bootstrap_threading(nsteps);
    bs = stats.BootstrapAverage();
    }
    //Write output
    WriteOutputs write;
    write.write(outfile,probfile,bw,bs,nsteps);
	return 0;
};


