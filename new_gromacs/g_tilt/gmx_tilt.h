#include <iostream>
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


typedef struct
{
    gmx_ana_selection_t *refsel;
    FILE                *fp;
    real                *ave;
    real                *n;
    gmx_ana_nbsearch_t  *nb;
} t_analysisdata;

int gmx_tilt(int argc, char *argv[]);
static int analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc, 
                         int nr, gmx_ana_selection_t *sel[], void *data);

