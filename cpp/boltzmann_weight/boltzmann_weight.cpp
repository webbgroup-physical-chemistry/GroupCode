#include "src/MyOptionParser.h"


using namespace std; 

int main(int argc, char * argv[])
{
	OptionParser parser;
	parser.assign_header(\
	"This is intended to calculate a Boltzmann weigghted average\n\
	given some value over instances, with each instance assigned\n\
	either a bin or a probability.  If each instance is assigned\n\
	a bin, a file governing the probability of each bin is also\n\
	required");
	parser.read_cmdline(argv, argv+argc);

	return 0;
};


