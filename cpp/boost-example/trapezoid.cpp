#include "integration.h"
/* Integration using the trapezoid rule */

double trapezoid(double (*f)(double),double x, int nbins)
{
    double dx = x / nbins;
    double sum = 0.;
    double sumsq = 0.;
    for (int i=0; i<nbins; i++){
        sum += (f(i*dx)+f(i*dx+dx))/2*dx;
    }
    std::cout << "\n Trapezoidal Integration using " << nbins << " bins." << std::endl;
    std::cout << "\t f(x) = " << sum << std::endl;
    return 0;
}
