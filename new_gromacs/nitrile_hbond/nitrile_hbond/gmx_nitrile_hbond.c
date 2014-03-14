//
//  gmx_nitrile_hbond.c
//  nitrile_hbond
//
//  Created by Andrew Ritchie on 3/5/14.
//  Copyright (c) 2014 Andrew Ritchie. All rights reserved.
//

#include "gmx_nitrile_hbond.h"


void gmx_nitrile_hbond(int argc, char *argv[])
{
    const char          *desc[] = {
        "Find waters hydrogen bonding to a specific atom.",
    };
    // Command line args
    int                 atom;
    gmx_bool            bVerbose;
    t_pargs             pa[] = {
        {"-atom", TRUE, etINT, {&atom}},
        {"-v", FALSE, etBOOL, {&bVerbose}},
    };
    // Output file
    t_filenm            fnm[] = {
        {efXVG, "-o", "hbonds", ffOPTWR},
    };
    
#define NFILE asize(fnm)
    
    gmx_ana_traj_t      *trj;
    output_env_t        oenv;
    t_analysisdata      d;
    int                 ngrps;
    gmx_ana_selection_t **sel;
    int                 g;
    int                 rc;
    
    CopyRight(stderr,argv[0]);
    /* Here, we can use flags to specify requirements for the selections and/or
     * other features of the library. */
    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nrefgrps(trj, 1);
    gmx_ana_set_nanagrps(trj, -1);
    
    /* If required, other functions can also be used to configure the library
     * before calling parse_trjana_args(). */
    
    // it's seg faulting here
    parse_trjana_args(trj, &argc, argv, PCA_CAN_VIEW,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);
    
    /* You can now do any initialization you wish, using the information
     * from the trj structure.
     * In particular, you should store any command-line parameter values that
     * the analysis part requires into d for them to be accessible. */
    
    /* First, we get some selection information from the structure */
    gmx_ana_get_refsel(trj, 0, &d.refsel);
    gmx_ana_get_nanagrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);
    
    /* First, we initialize the neighborhood search for the first index
     * group. */
    return;
}