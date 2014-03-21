#include <iostream>
#include <vector>
#include <cmath>
// gromacs c headers
#ifdef __cplusplus
extern "C"  
#endif
#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/nbsearch.h>
#include <gromacs/trajana.h>
#include <gromacs/tpxio.h>


typedef struct
{
    FILE   *fp;
    std::vector<double> distance;
    std::vector<double> x_rotation;
    std::vector<double> y_rotation;
    std::vector<double> z_rotation;
    std::vector<double> theta_rotation;
    std::vector<double> phi_rotation;
} t_tiltdata;

typedef std::vector< std::vector<double> > Matrix;

int gmx_tilt(int argc, char *argv[]);
void readCoords(int n_atoms, atom_id ind_atoms[], t_trxframe *fr, t_topology *top, 
                Matrix &coords, gmx_bool bVerbose);
void readCoords(int n_atoms, atom_id ind_atoms[], rvec *x, t_topology *top, 
                   Matrix &coords, gmx_bool bVerbose);
void displacement( Matrix &reference, Matrix &frame, t_tiltdata  &data );
void average_coordinate( Matrix &coords, std::vector<double> &xyz);
