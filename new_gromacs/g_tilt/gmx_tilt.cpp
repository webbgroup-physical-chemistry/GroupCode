#include "gmx_tilt.h"

int gmx_tilt(int argc, char *argv[])
{
    const char      *desc[] = {
        "Look at the tilt angle of RalGDS when docked to different",
        "GTPases: p21Ras and Rap1a",
    };

    gmx_bool        bVerbose = FALSE;
    const char      *struct_file, *traj_file;
    const char      *ndx_file, *xvg_file;
    const char      *ref_file;
    t_tiltdata  data;
    t_pargs         pa[] = {
        { "-v", FALSE, etBOOL, {&bVerbose},
          "Be slightly more verbose"},
    };
    t_filenm        fnm[] = {
        {efTPS, NULL, NULL, ffREAD},
        {efTRX, NULL, NULL, ffREAD},
        //{efGRO, "-fr", NULL, ffREAD},
        {efTPS, "-fr", NULL, ffREAD},
        {efXVG, "-o","tilt",ffOPTWR},
        {efNDX, NULL, NULL, ffREAD},
    };
#define NFILE asize(fnm)
#define NPA asize(pa)
    output_env_t    oenv;
    int             ngrps, nrefgrps;
    t_topology      struct_top, ref_top;
    t_atoms         *struct_atoms=NULL, *ref_atoms=NULL;
    t_trxframe      struct_fr, ref_fr;
    t_trxstatus     *struct_status, *ref_status;
    rvec            *struct_xtop, *ref_xtop;
    matrix          struct_box, ref_box;
    int             struct_ePBC, ref_ePBC;
    int             struct_flags=TRX_READ_X, ref_flags=TRX_READ_X;
    char            buffer[1024];

    CopyRight(stderr,argv[0]);
    parse_common_args(&argc, argv, 
                      PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW | 
                      PCA_TIME_UNIT | PCA_BE_NICE | PCA_CAN_TIME,
                      NFILE, fnm, NPA, pa, asize(desc), desc,
                      0, NULL, &oenv);

    /* Get inputs */
    struct_file = ftp2fn(efTPS, NFILE, fnm);
    traj_file   = opt2fn( "-f", NFILE, fnm);
    ref_file    = opt2fn("-fr", NFILE, fnm); 
    ndx_file    = ftp2fn(efNDX, NFILE, fnm);
    
    /* Open inputs */
    read_tps_conf(struct_file, buffer, &struct_top, &struct_ePBC,
                  &struct_xtop, NULL, struct_box, TRUE);
    sfree(struct_xtop);
    struct_atoms = &struct_top.atoms;

    read_tps_conf(ref_file, buffer, &ref_top, &ref_ePBC,
                  &ref_xtop, NULL, ref_box, TRUE);
    ref_atoms = &ref_top.atoms;

    /* open xvg file */
    data.fp = NULL;
    if (opt2bSet("-o", NFILE, fnm))
    {
        data.fp = xvgropen(opt2fn("-o", NFILE, fnm), "Center of Mass displacement","Time[ps]", "Distance[nm]", oenv);
    }
    
    /* Get index information */
    static int      n_structatoms, n_refatoms;
    atom_id         *ind_structatoms, *ind_refatoms;
    char            *gn_structatoms, *gn_refatoms;
    int             i_structatoms, i_refatoms;
    std::cout<< "Select group from <" << ref_file << "> to compare to." << std::endl;
    get_index(ref_atoms, ndx_file, 1, &n_refatoms, &ind_refatoms, &gn_refatoms);
    std::cout<< "Select group from <" << struct_file <<"> to compare to <" << ref_file << ">." << std::endl;
    get_index(struct_atoms, ndx_file, 1, &n_structatoms, &ind_structatoms, &gn_structatoms);

    if ( n_structatoms != n_refatoms ) 
    {
        gmx_fatal(FARGS, "\nThe number of atoms in the index groups are different.");
    }
     
    /* read the reference file and get coordinates */
    Matrix ref_ndxatoms(n_refatoms);
    readCoords(n_refatoms, ind_refatoms, ref_xtop, &ref_top, ref_ndxatoms, bVerbose);

    /* read trajectory file and loop through frames */
    int i=0;
    read_first_frame(oenv, &struct_status, traj_file, &struct_fr, struct_flags);
    do {
        Matrix struct_ndxatoms(n_structatoms);
        readCoords(n_structatoms, ind_structatoms, &struct_fr, &struct_top, struct_ndxatoms, bVerbose);
        displacement(ref_ndxatoms, struct_ndxatoms, data);
        if (bVerbose)
        {
            std::cout <<"\n Distance: " << data.distance[i];
        }
/*
1. Write kabsch algorythm
2. Figure out how to convert rotation matrix into angle
3. Convert angles to spherical coordinates
*/

        if (bVerbose)
        {
            std::cout << "\n";
        }
        i++;
    } while(read_next_frame(oenv, struct_status, &struct_fr));



    return 0;
}

void average_coordinate( Matrix &coords, std::vector<double> &xyz)
{
    for (int i=0; i<coords.size(); i++)
    {
        for (int j=0; j<coords[0].size(); j++)
        {
            xyz[j] += coords[i][j];
        }
    }
    return;
}

void displacement( Matrix &reference, Matrix &frame, t_tiltdata &data)
{
    std::vector<double> reference_xyz(3,0), frame_xyz(3,0);
    average_coordinate(reference, reference_xyz);
    average_coordinate(frame, frame_xyz);
    double x, y, z;
    x = reference_xyz[0] - frame_xyz[0];
    y = reference_xyz[1] - frame_xyz[1];
    z = reference_xyz[2] - frame_xyz[2];
    data.distance.push_back( sqrt(x*x+y*y+z*z) );
}

/* Reading a TPX */
void readCoords(int n_atoms, atom_id ind_atoms[], rvec *x, t_topology *top, Matrix &coords, gmx_bool bVerbose)
{
    double mass = 0;
    for (int i=0; i<n_atoms; i++)
    {
        int n = ind_atoms[i];
        std::vector<double> xyz(3);
        xyz[0] = x[n][0];
        xyz[1] = x[n][1];
        xyz[2] = x[n][2];
        mass += top->atoms.atom[n].m;
        coords[i] = xyz;
    }
    for (int i=0; i<n_atoms; i++)
    {
        int n = ind_atoms[i];
        for (int j=0; j<3; j++)
        {
            coords[i][j] *= top->atoms.atom[n].m;
            coords[i][j] /= mass;
        }
    }
    if (bVerbose)
    {
        std::cout << "\n\n Mass Weighted Atom Coordinates: ";
        for (int i=0; i<n_atoms; i++)
        {
            std::cout << "\n\t" << *top->atoms.atomname[ind_atoms[i]];
            std::cout << "\t" << top->atoms.atom[ind_atoms[i]].m << " " << coords[i][0] << " " << coords[i][1] << " " << coords[i][2];
       }
        std::cout << "\n";
    }
    return;
}
/* Reading a trajectory */
void readCoords(int n_atoms, atom_id ind_atoms[], t_trxframe *fr, t_topology *top, Matrix &coords, gmx_bool bVerbose)
{
    double mass = 0;
    for (int i=0; i<n_atoms; i++)
    {
        int n = ind_atoms[i];
        std::vector<double> xyz(3);
        xyz[0] = fr->x[n][XX];
        xyz[1] = fr->x[n][YY];
        xyz[2] = fr->x[n][ZZ];
        mass += top->atoms.atom[n].m;
        coords[i] = xyz;
    }
    for (int i=0; i<n_atoms; i++)
    {
        int n = ind_atoms[i];
        for (int j=0; j<3; j++)
        {
            coords[i][j] *= top->atoms.atom[n].m;
            coords[i][j] /= mass;
        }
    }
    if (bVerbose)
    {
        std::cout << "\n\n Mass Weighted Atom Coordinates: ";
        for (int i=0; i<n_atoms; i++)
        {
            std::cout << "\n\t" << *top->atoms.atomname[ind_atoms[i]];
            std::cout << "\t" << top->atoms.atom[ind_atoms[i]].m << " " << coords[i][0] << " " << coords[i][1] << " " << coords[i][2];
       }
        std::cout << "\n";
    }
    return;
}

