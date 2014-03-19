#include "gmx_tilt.h"

int gmx_tilt(int argc, char *argv[])
{
    const char      *desc[] = {
        "Look at the tilt angle of RalGDS when docked to different",
        "GTPases: p21Ras and Rap1a",
    };

    gmx_bool        bVerbose = FALSE;
    const char      *struct_file, *traj_file, *ndx_file, *xvg_file;
    t_pargs         pa[] = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be slightly more verbose"},
    };
    t_filenm        fnm[] = {
        {efTPS, NULL, NULL, ffREAD},
        {efTRX, NULL, NULL, ffREAD},
        {efTPS, "-fr", NULL, ffREAD},
        {efXVG, "-o","tilt",ffOPTWR},
        {efNDX, NULL, NULL, ffREAD},
    };
#define NFILE asize(fnm)
#define NPA asize(pa)
    output_env_t    oenv;
    int             ngrps, nrefgrps;

    CopyRight(stderr,argv[0]);
    parse_common_args(&argc, argv, 
                      PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW | 
                      PCA_TIME_UNIT | PCA_BE_NICE | PCA_CAN_TIME,
                      NFILE, fnm, NPA, pa, asize(desc), desc,
                      0, NULL, &oenv);

    return 0;
}


/*! \brief
 * Function that does the analysis for a single frame.
 *
 * It is called once for each frame.
 */
static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata     *d = (t_analysisdata *)data;
    int                 g, i;
    real                frave;
    int                 rc;

    /* Here, you can do whatever analysis your program requires for a frame. */
    if (d->fp)
    {
        fprintf(d->fp, "%10.3f", fr->time);
    }
    rc = gmx_ana_nbsearch_pos_init(d->nb, pbc, &d->refsel->p);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "Neighborhood search initialization failed");
    }
    for (g = 0; g < nr; ++g)
    {
        frave = 0;
        for (i = 0; i < sel[g]->p.nr; ++i)
        {
            frave += gmx_ana_nbsearch_pos_mindist(d->nb, &sel[g]->p, i);
        }
        d->ave[g] += frave;
        d->n[g]   += sel[g]->p.nr;
        if (d->fp)
        {
            frave /= sel[g]->p.nr;
            fprintf(d->fp, " %.3f", frave);
        }
    }
    if (d->fp)
    {
        fprintf(d->fp, "\n");
    }
    /* We need to return 0 to tell that everything went OK */
    return 0;
}

/*! \brief
 * Function that implements the analysis tool.
 *
 * Following the style of Gromacs analysis tools, this function is called
 * \p gmx_something.
 */
