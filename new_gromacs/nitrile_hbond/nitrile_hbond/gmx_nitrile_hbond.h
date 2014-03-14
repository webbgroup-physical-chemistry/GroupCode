//
//  gmx_nitrile_hbond.h
//  nitrile_hbond
//
//  Created by Andrew Ritchie on 3/5/14.
//  Copyright (c) 2014 Andrew Ritchie. All rights reserved.
//

#ifndef nitrile_hbond_gmx_nitrile_hbond_h
#define nitrile_hbond_gmx_nitrile_hbond_h

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

#endif

typedef struct
{
    gmx_ana_selection_t *refsel;
    FILE                *fp;
    real                *ave;
    real                *n;
    gmx_ana_nbsearch_t  *nb;
} t_analysisdata;

void gmx_nitrile_hbond(int argc, char *argv[]);
