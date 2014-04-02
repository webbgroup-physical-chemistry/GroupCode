#include "integration.h"
/* Integration using the trapezoid rule */

double simpson(double (*f)(double),double x, int nbins)
{
    double dx = x / nbins;
    double sum = 0.;
    double sumsq = 0.;
    double a, b;
    for (int i=0; i<nbins; i++){
        a = i * dx;
        b = a + dx;
        sum += (b-a)/6*(f(a)+4*f((a+b)/2)+f(b));
    }

    std::cout << "\n Simpson's Rule Integration using " << nbins << " bins." << std::endl;
    std::cout << "\t f(x) = " << sum << std::endl;
    return 0;
}
