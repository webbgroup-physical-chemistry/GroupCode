#include "integration.h"
#include <boost/program_options.hpp>
#include <iterator>

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    double tbins,sbins;
    double mctol;
    try {
        po::options_description desc("Options");
        desc.add_options()
            ("help,h", "This stuff!")
            ("trapezoid_bins,t", po::value<double>(&tbins)->default_value(1e2,"100"),
                "Number of bins for trapezoidal integration")
            ("simpon_bins,s", po::value<double>(&sbins)->default_value(1e2,"100"),
                "Number of bins for Simopson's rule integration")
            ("mc_tol,m", po::value<double>(&mctol)->default_value(1e-3),
                "Monte Carlo integration tolerance")
        ;
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
    }
    catch(...) {
       std::cerr << "Exception of unknown type!\n";
    }

    trapezoid(func,1,tbins);
    simpson(func,1,sbins);
    mc(func,1,mctol);
    return 0;
}
