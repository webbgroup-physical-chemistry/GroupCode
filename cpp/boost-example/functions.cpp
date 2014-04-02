#include "integration.h"

/* 
    erf(x) = 2/sqrt[pi] INT[ exp[-u^2] du, 0, x ] 
    so we want f(x) = 2/sqrt[pi] * exp[-x*x]    
*/

double func(double x)
{
    return 2 / sqrt(M_PI) * std::exp(-x*x);
}

double t_func(double x)
{
    return cos(x);
}


