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
                       
	parser.parse_args();
	string probfile = parser.Options()[0];
	string listfile = parser.Options()[1];

	Probability probability;
	probability.ReadProb( probfile );
	std::vector<double> prob = probability.Prob();

	std::cout << "There are " << prob.size()-1 << " bins." << std::endl;
/*	double total = 0;
	for (int i = 0 ; i < prob.size() ; i++ ) 
	{
  		std::cout << i << " " << prob[i] << std::endl;
		total += prob[i];
	}
	std::cout << total << std::endl;
	std::cout << prob.size() << std::endl;
*/

//	Dats b;
//	b.ReadDat( blah );
//	std::vector<int> bb = b.Bin();
//	for (int i=0;i<bb.size();i++){std::cout<<i<<" "<<bb[i]<<"\n";}

	return 0;
};


