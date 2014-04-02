#include "integration.h"
/* Integration using monte carlo method */

double gq2(double (*f)(double), double x, double tol)
//double mc(double x, double tol)
{
    /* set up random number generator */
    const gsl_rng_type *T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    double seed = time(NULL);
    gsl_rng_set(r,seed);
    /* Int[f(x)dx,a,b] = (b-a)*f_average */
    /* MC part */
    if (tol > 1){ tol = 1e-3; }
    double u, fu;
    double result;
    double error = 1.;
    int nstep = 0;
    double sum = 0.0, sumsq = 0.0;

    while (error > tol) {
        u = gsl_rng_uniform(r)*x; // 0 + random * (high - low)
        fu = f(u);
        sum += fu;
        sumsq += fu*fu;
        nstep ++;

        result = sum/nstep * (x-0);
        if (nstep > 10){
            error = sqrt((sumsq/nstep - sum*sum/nstep/nstep)/nstep) * (x-0);
        }
    }
    result = floorf(result/tol)*tol;
    std::cout << "\nConverged after " << nstep << " steps:\n " << result << " +/- " << error << std::endl;
    std::cout << "\n Range: \t" << result-error << " <= f(x) <= " << result+error << std::endl;

    return 0;
}

