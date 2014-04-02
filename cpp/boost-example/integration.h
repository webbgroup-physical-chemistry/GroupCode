#include <iostream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI = 4*atan(1)
#endif


double func(double x);
double t_func(double x);
double mc(double (*f)(double),double x,double tol);
double trapezoid(double (*f)(double),double x,int nbins);
double simpson(double (*f)(double), double x,int nbins);
double gq();


